################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../apps/GSPMHD2D/init.cpp \
../apps/GSPMHD2D/integral.cpp \
../apps/GSPMHD2D/main.cpp \
../apps/GSPMHD2D/util.cpp 

O_SRCS += \
../apps/GSPMHD2D/init.o \
../apps/GSPMHD2D/integral.o \
../apps/GSPMHD2D/main.o \
../apps/GSPMHD2D/util.o 

OBJS += \
./apps/GSPMHD2D/init.o \
./apps/GSPMHD2D/integral.o \
./apps/GSPMHD2D/main.o \
./apps/GSPMHD2D/util.o 

CPP_DEPS += \
./apps/GSPMHD2D/init.d \
./apps/GSPMHD2D/integral.d \
./apps/GSPMHD2D/main.d \
./apps/GSPMHD2D/util.d 


# Each subdirectory must supply rules for building sources it contributes
apps/GSPMHD2D/%.o: ../apps/GSPMHD2D/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/guo/Research/SimulationCode/FDPS-master/src" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


