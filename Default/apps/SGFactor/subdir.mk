################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../apps/SGFactor/SGFactorMesh.cpp 

O_SRCS += \
../apps/SGFactor/SGFactorMesh.o 

OBJS += \
./apps/SGFactor/SGFactorMesh.o 

CPP_DEPS += \
./apps/SGFactor/SGFactorMesh.d 


# Each subdirectory must supply rules for building sources it contributes
apps/SGFactor/%.o: ../apps/SGFactor/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/guo/Research/SimulationCode/FDPS-master/src" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


