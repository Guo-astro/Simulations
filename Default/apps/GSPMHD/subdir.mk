################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../apps/GSPMHD/init.cpp \
../apps/GSPMHD/integral.cpp \
../apps/GSPMHD/main.cpp \
../apps/GSPMHD/util.cpp 

O_SRCS += \
../apps/GSPMHD/init.o \
../apps/GSPMHD/integral.o \
../apps/GSPMHD/main.o \
../apps/GSPMHD/util.o 

OBJS += \
./apps/GSPMHD/init.o \
./apps/GSPMHD/integral.o \
./apps/GSPMHD/main.o \
./apps/GSPMHD/util.o 

CPP_DEPS += \
./apps/GSPMHD/init.d \
./apps/GSPMHD/integral.d \
./apps/GSPMHD/main.d \
./apps/GSPMHD/util.d 


# Each subdirectory must supply rules for building sources it contributes
apps/GSPMHD/%.o: ../apps/GSPMHD/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/guo/Research/SimulationCode/FDPS-master/src" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


