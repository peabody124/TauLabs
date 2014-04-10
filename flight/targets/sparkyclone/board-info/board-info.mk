BOARD_TYPE          := 0x88
BOARD_REVISION      := 0x02
BOOTLOADER_VERSION  := 0x81
HW_TYPE             := 0x00		# seems to be unused

MCU                 := cortex-m4
CHIP                := STM32F303CCT
BOARD               := STM32F30x_SPARKY
MODEL               := HD
MODEL_SUFFIX        := 

OPENOCD_JTAG_CONFIG := stlink-v2.cfg
OPENOCD_CONFIG      := stm32f3xx.stlink.cfg

FW_BANK_BASE        := 0x08000000  # Start of firmware flash @48kB
FW_BANK_SIZE        := 0x00038000  # Should include FW_DESC_SIZE (208kB)

FW_DESC_SIZE        := 0x00000064

EE_BANK_BASE        := 0x08038000
EE_BANK_SIZE        := 0x00008000  # (32kb)

EF_BANK_BASE        := 0x08000000  # Start of entire flash image (usually start of bootloader as well)
EF_BANK_SIZE        := 0x00040000  # Size of the entire flash image (from bootloader until end of firmware)

OSCILLATOR_FREQ     :=   8000000
SYSCLK_FREQ         :=  72000000

