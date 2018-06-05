################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../apps/nbodysph/init.cpp \
../apps/nbodysph/integral.cpp \
../apps/nbodysph/main.cpp \
../apps/nbodysph/util.cpp 

OBJS += \
./apps/nbodysph/init.o \
./apps/nbodysph/integral.o \
./apps/nbodysph/main.o \
./apps/nbodysph/util.o 

CPP_DEPS += \
./apps/nbodysph/init.d \
./apps/nbodysph/integral.d \
./apps/nbodysph/main.d \
./apps/nbodysph/util.d 


# Each subdirectory must supply rules for building sources it contributes
apps/nbodysph/%.o: ../apps/nbodysph/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/guo/Research/SimulationCode/FDPS-master/src" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


