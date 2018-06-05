################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tests/DomainInfo/decomposeDomain/mainf32.cpp \
../tests/DomainInfo/decomposeDomain/mainf64.cpp 

OBJS += \
./tests/DomainInfo/decomposeDomain/mainf32.o \
./tests/DomainInfo/decomposeDomain/mainf64.o 

CPP_DEPS += \
./tests/DomainInfo/decomposeDomain/mainf32.d \
./tests/DomainInfo/decomposeDomain/mainf64.d 


# Each subdirectory must supply rules for building sources it contributes
tests/DomainInfo/decomposeDomain/%.o: ../tests/DomainInfo/decomposeDomain/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/guo/Research/SimulationCode/FDPS-master/src" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


