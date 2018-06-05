################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/phantom_grape_x86/G5/newton/cpu.c \
../src/phantom_grape_x86/G5/newton/direct.c 

OBJS += \
./src/phantom_grape_x86/G5/newton/cpu.o \
./src/phantom_grape_x86/G5/newton/direct.o 

C_DEPS += \
./src/phantom_grape_x86/G5/newton/cpu.d \
./src/phantom_grape_x86/G5/newton/direct.d 


# Each subdirectory must supply rules for building sources it contributes
src/phantom_grape_x86/G5/newton/%.o: ../src/phantom_grape_x86/G5/newton/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


