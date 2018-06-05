################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/phantom_grape_x86/G5/table/gravity_kernel.c \
../src/phantom_grape_x86/G5/table/gravity_kernel_avx2.c \
../src/phantom_grape_x86/G5/table/gravity_kernel_avx2_ntdr.c \
../src/phantom_grape_x86/G5/table/pg5_table.c \
../src/phantom_grape_x86/G5/table/phantom_g5.c 

OBJS += \
./src/phantom_grape_x86/G5/table/gravity_kernel.o \
./src/phantom_grape_x86/G5/table/gravity_kernel_avx2.o \
./src/phantom_grape_x86/G5/table/gravity_kernel_avx2_ntdr.o \
./src/phantom_grape_x86/G5/table/pg5_table.o \
./src/phantom_grape_x86/G5/table/phantom_g5.o 

C_DEPS += \
./src/phantom_grape_x86/G5/table/gravity_kernel.d \
./src/phantom_grape_x86/G5/table/gravity_kernel_avx2.d \
./src/phantom_grape_x86/G5/table/gravity_kernel_avx2_ntdr.d \
./src/phantom_grape_x86/G5/table/pg5_table.d \
./src/phantom_grape_x86/G5/table/phantom_g5.d 


# Each subdirectory must supply rules for building sources it contributes
src/phantom_grape_x86/G5/table/%.o: ../src/phantom_grape_x86/G5/table/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


