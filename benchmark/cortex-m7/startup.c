/*
Bare-metal support for the NMPC Cortex-M7 benchmark under QEMU mps2-an500. Provides the vector
table, start-up, semihosting output, an instruction counter built on SysTick
(QEMU -icount shift=0 advances virtual time by 1 ns per instruction, and the
25 MHz SysTick therefore ticks once per 40 instructions), and a stack
high-water mark.
*/
#include <stdint.h>

extern uint32_t _estack, _sstack, _sidata, _sdata, _edata, _sbss, _ebss;
int main(void);
void __libc_init_array(void);

#define SYST_CSR (*(volatile uint32_t*)0xE000E010)
#define SYST_RVR (*(volatile uint32_t*)0xE000E014)
#define SYST_CVR (*(volatile uint32_t*)0xE000E018)
#define CPACR    (*(volatile uint32_t*)0xE000ED88)

static volatile uint32_t wraps;

static int semihost(int op, void* arg) {
    register int r0 __asm__("r0") = op;
    register void* r1 __asm__("r1") = arg;
    __asm__ volatile("bkpt 0xAB" : "+r"(r0) : "r"(r1) : "memory");
    return r0;
}

void print(const char* s) { semihost(0x04, (void*)s); }

void print_u(uint64_t v) {
    char buffer[24]; int i = 23; buffer[i] = 0;
    do { buffer[--i] = (char)('0' + v % 10); v /= 10; } while (v);
    print(&buffer[i]);
}

void finish(void) { semihost(0x18, (void*)0x20026); for (;;) {} }

void SysTick_Handler(void) { wraps++; }
void Fault_Handler(void) { print("FAULT\n"); finish(); }

void timer_start(void) { SYST_RVR = 0xFFFFFF; SYST_CVR = 0; SYST_CSR = 7; }

/* Instructions executed since timer_start (resolution 40). */
uint64_t instructions(void) {
    uint32_t before, after, value;
    do { before = wraps; value = SYST_CVR; after = wraps; } while (before != after);
    return (((uint64_t)before << 24) + (0xFFFFFFu - value)) * 40u;
}

static uint32_t* paint_top;

__attribute__((noinline)) void stack_paint(void) {
    uint32_t* top = (uint32_t*)__builtin_frame_address(0) - 32;
    for (uint32_t* p = &_sstack; p < top; p++) *p = 0xA5A5A5A5u;
    paint_top = top;
}

/*
Bytes of stack in use at the deepest point since stack_paint was called, counted from the stack's
original top (the reset value of the stack pointer), so that it includes start-up's and main's own
frames as well as everything called below them.
*/
uint32_t stack_used(void) {
    uint32_t* p = &_sstack;
    while (p < paint_top && *p == 0xA5A5A5A5u) p++;
    return (uint32_t)((uint8_t*)&_estack - (uint8_t*)p);
}

void _init(void) {}

void Reset_Handler(void) {
    CPACR |= (0xFu << 20);      /* enable the FPU */
    __asm__ volatile("dsb\n isb");
    uint32_t* source = &_sidata;
    for (uint32_t* p = &_sdata; p < &_edata; p++) *p = *source++;
    for (uint32_t* p = &_sbss; p < &_ebss; p++) *p = 0;
    __libc_init_array();
    main();
    finish();
}

__attribute__((section(".vectors"), used))
void (* const vectors[16])(void) = {
    (void (*)(void))&_estack, Reset_Handler, Fault_Handler, Fault_Handler, Fault_Handler, Fault_Handler,
    Fault_Handler, 0, 0, 0, 0, Fault_Handler, Fault_Handler, 0, Fault_Handler, SysTick_Handler
};
