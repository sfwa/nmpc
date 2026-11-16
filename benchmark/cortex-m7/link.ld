/* Linker script for the NMPC Cortex-M7 benchmark under QEMU mps2-an500 (Cortex-M7): 4 MB code at 0x00000000, 4 MB data at 0x20000000. */
ENTRY(Reset_Handler)
MEMORY {
    ROM (rx)  : ORIGIN = 0x00000000, LENGTH = 4M
    RAM (rwx) : ORIGIN = 0x20000000, LENGTH = 4M
}
_estack = ORIGIN(RAM) + LENGTH(RAM);
_sstack = _estack - 0x80000;    /* 512 kB reserved so that stack use can be measured */
SECTIONS {
    .text : {
        KEEP(*(.vectors))
        *(.text*) *(.rodata*)
        . = ALIGN(4);
        __preinit_array_start = .; KEEP(*(.preinit_array*)) __preinit_array_end = .;
        __init_array_start = .; KEEP(*(SORT(.init_array.*))) KEEP(*(.init_array*)) __init_array_end = .;
        __fini_array_start = .; KEEP(*(.fini_array*)) __fini_array_end = .;
    } > ROM
    .ARM.extab : { *(.ARM.extab* .gnu.linkonce.armextab.*) } > ROM
    .ARM.exidx : { __exidx_start = .; *(.ARM.exidx*) __exidx_end = .; } > ROM
    _sidata = LOADADDR(.data);
    .data : { _sdata = .; *(.data*) . = ALIGN(4); _edata = .; } > RAM AT > ROM
    .bss (NOLOAD) : { _sbss = .; *(.bss*) *(COMMON) . = ALIGN(4); _ebss = .; end = .; } > RAM
}
