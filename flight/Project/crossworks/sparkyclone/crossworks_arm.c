#include <pios.h>

//libc_user_libc_v7em_fpv4_sp_d16_hard_t_le_eabi.a(user_libc.o): In function `__aeabi_read_tp':
//user_libc.c:(.text.libc.__aeabi_read_tp+0x4): undefined reference to `__tbss_start__'
volatile void* __tbss_start__;

// ---------------------------------------------------------
#if 1

typedef struct _hardfault_args_t {
  unsigned int r0;
  unsigned int r1;
  unsigned int r2;
  unsigned int r3;
  unsigned int r12;
  unsigned int lr;
  unsigned int pc;
  unsigned int psr;
} hardfault_args_t;

void hard_fault_handler_c (hardfault_args_t* hardfault_args)
{
   for(;;) {
     __asm("BKPT #0\n") ;
   }
}

void HardFault_Handler(void) __attribute__((naked));
void HardFault_Handler(void)
{
    __asm("TST LR, #4\r\n"
        "ITE EQ\r\n"
        "MRSEQ R0, MSP\r\n"
        "MRSNE R0, PSP\r\n"
        "B hard_fault_handler_c");
}

void NMI_Handler(void) __attribute__((naked));
void NMI_Handler(void)
{

  for(;;) {
    __asm("BKPT #0\n") ;
  }
}

void MemManage_Handler(void) __attribute__((naked));
void MemManage_Handler(void)
{

  for(;;) {
    __asm("BKPT #0\n") ;
  }
}

void BusFault_Handler(void) __attribute__((naked));
void BusFault_Handler(void)
{

  for(;;) {
    __asm("BKPT #0\n") ;
  }
}


void UsageFault_Handler(void) __attribute__((naked));
void UsageFault_Handler(void)
{

  for(;;) {
    __asm("BKPT #0\n") ;
  }
}

void WWDG_IRQHandler(void) __attribute__((naked));
void WWDG_IRQHandler(void)
{

  for(;;) {
    __asm("BKPT #0\n") ;
  }
}

#endif
