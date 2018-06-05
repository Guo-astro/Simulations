################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../apps/HVCC_CHEM_MBH5_EXP/init.cpp \
../apps/HVCC_CHEM_MBH5_EXP/integral.cpp \
../apps/HVCC_CHEM_MBH5_EXP/main.cpp \
../apps/HVCC_CHEM_MBH5_EXP/util.cpp 

O_SRCS += \
../apps/HVCC_CHEM_MBH5_EXP/integral.o \
../apps/HVCC_CHEM_MBH5_EXP/main.o \
../apps/HVCC_CHEM_MBH5_EXP/util.o 

OBJS += \
./apps/HVCC_CHEM_MBH5_EXP/init.o \
./apps/HVCC_CHEM_MBH5_EXP/integral.o \
./apps/HVCC_CHEM_MBH5_EXP/main.o \
./apps/HVCC_CHEM_MBH5_EXP/util.o 

CPP_DEPS += \
./apps/HVCC_CHEM_MBH5_EXP/init.d \
./apps/HVCC_CHEM_MBH5_EXP/integral.d \
./apps/HVCC_CHEM_MBH5_EXP/main.d \
./apps/HVCC_CHEM_MBH5_EXP/util.d 


# Each subdirectory must supply rules for building sources it contributes
apps/HVCC_CHEM_MBH5_EXP/%.o: ../apps/HVCC_CHEM_MBH5_EXP/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/guo/Research/SimulationCode/FDPS-master/src" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


