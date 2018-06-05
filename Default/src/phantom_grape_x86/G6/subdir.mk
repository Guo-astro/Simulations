################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/phantom_grape_x86/G6/banana01.c \
../src/phantom_grape_x86/G6/calc_energy.c \
../src/phantom_grape_x86/G6/hermite.c \
../src/phantom_grape_x86/G6/io.c \
../src/phantom_grape_x86/G6/timeprof.c \
../src/phantom_grape_x86/G6/timestep.c 

OBJS += \
./src/phantom_grape_x86/G6/banana01.o \
./src/phantom_grape_x86/G6/calc_energy.o \
./src/phantom_grape_x86/G6/hermite.o \
./src/phantom_grape_x86/G6/io.o \
./src/phantom_grape_x86/G6/timeprof.o \
./src/phantom_grape_x86/G6/timestep.o 

C_DEPS += \
./src/phantom_grape_x86/G6/banana01.d \
./src/phantom_grape_x86/G6/calc_energy.d \
./src/phantom_grape_x86/G6/hermite.d \
./src/phantom_grape_x86/G6/io.d \
./src/phantom_grape_x86/G6/timeprof.d \
./src/phantom_grape_x86/G6/timestep.d 


# Each subdirectory must supply rules for building sources it contributes
src/phantom_grape_x86/G6/%.o: ../src/phantom_grape_x86/G6/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


