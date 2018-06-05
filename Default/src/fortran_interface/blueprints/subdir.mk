################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/fortran_interface/blueprints/FDPS_Manipulators_blueprint.cpp \
../src/fortran_interface/blueprints/FDPS_ftn_if_blueprint.cpp \
../src/fortran_interface/blueprints/main_blueprint.cpp 

OBJS += \
./src/fortran_interface/blueprints/FDPS_Manipulators_blueprint.o \
./src/fortran_interface/blueprints/FDPS_ftn_if_blueprint.o \
./src/fortran_interface/blueprints/main_blueprint.o 

CPP_DEPS += \
./src/fortran_interface/blueprints/FDPS_Manipulators_blueprint.d \
./src/fortran_interface/blueprints/FDPS_ftn_if_blueprint.d \
./src/fortran_interface/blueprints/main_blueprint.d 


# Each subdirectory must supply rules for building sources it contributes
src/fortran_interface/blueprints/%.o: ../src/fortran_interface/blueprints/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/guo/Research/SimulationCode/FDPS-master/src" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


