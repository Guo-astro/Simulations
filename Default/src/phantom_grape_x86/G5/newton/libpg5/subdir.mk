################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/phantom_grape_x86/G5/newton/libpg5/gravity.c \
../src/phantom_grape_x86/G5/newton/libpg5/gravity_avx2.c \
../src/phantom_grape_x86/G5/newton/libpg5/pg5_fortran.c \
../src/phantom_grape_x86/G5/newton/libpg5/phantom_g5.c \
../src/phantom_grape_x86/G5/newton/libpg5/rsqrt.c 

OBJS += \
./src/phantom_grape_x86/G5/newton/libpg5/gravity.o \
./src/phantom_grape_x86/G5/newton/libpg5/gravity_avx2.o \
./src/phantom_grape_x86/G5/newton/libpg5/pg5_fortran.o \
./src/phantom_grape_x86/G5/newton/libpg5/phantom_g5.o \
./src/phantom_grape_x86/G5/newton/libpg5/rsqrt.o 

C_DEPS += \
./src/phantom_grape_x86/G5/newton/libpg5/gravity.d \
./src/phantom_grape_x86/G5/newton/libpg5/gravity_avx2.d \
./src/phantom_grape_x86/G5/newton/libpg5/pg5_fortran.d \
./src/phantom_grape_x86/G5/newton/libpg5/phantom_g5.d \
./src/phantom_grape_x86/G5/newton/libpg5/rsqrt.d 


# Each subdirectory must supply rules for building sources it contributes
src/phantom_grape_x86/G5/newton/libpg5/%.o: ../src/phantom_grape_x86/G5/newton/libpg5/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


