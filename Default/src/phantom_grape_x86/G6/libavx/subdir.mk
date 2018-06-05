################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/phantom_grape_x86/G6/libavx/gravity.c \
../src/phantom_grape_x86/G6/libavx/gravity_avx2.c \
../src/phantom_grape_x86/G6/libavx/phantom_g6.c \
../src/phantom_grape_x86/G6/libavx/timeprof.c 

OBJS += \
./src/phantom_grape_x86/G6/libavx/gravity.o \
./src/phantom_grape_x86/G6/libavx/gravity_avx2.o \
./src/phantom_grape_x86/G6/libavx/phantom_g6.o \
./src/phantom_grape_x86/G6/libavx/timeprof.o 

C_DEPS += \
./src/phantom_grape_x86/G6/libavx/gravity.d \
./src/phantom_grape_x86/G6/libavx/gravity_avx2.d \
./src/phantom_grape_x86/G6/libavx/phantom_g6.d \
./src/phantom_grape_x86/G6/libavx/timeprof.d 


# Each subdirectory must supply rules for building sources it contributes
src/phantom_grape_x86/G6/libavx/%.o: ../src/phantom_grape_x86/G6/libavx/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


