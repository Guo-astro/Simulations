################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/particle_mesh/decomposition.cpp \
../src/particle_mesh/evolve.cpp \
../src/particle_mesh/io.cpp \
../src/particle_mesh/pm_parallel.cpp \
../src/particle_mesh/pp2.cpp \
../src/particle_mesh/ptobmp.cpp \
../src/particle_mesh/tools.cpp 

OBJS += \
./src/particle_mesh/decomposition.o \
./src/particle_mesh/evolve.o \
./src/particle_mesh/io.o \
./src/particle_mesh/pm_parallel.o \
./src/particle_mesh/pp2.o \
./src/particle_mesh/ptobmp.o \
./src/particle_mesh/tools.o 

CPP_DEPS += \
./src/particle_mesh/decomposition.d \
./src/particle_mesh/evolve.d \
./src/particle_mesh/io.d \
./src/particle_mesh/pm_parallel.d \
./src/particle_mesh/pp2.d \
./src/particle_mesh/ptobmp.d \
./src/particle_mesh/tools.d 


# Each subdirectory must supply rules for building sources it contributes
src/particle_mesh/%.o: ../src/particle_mesh/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"/Users/guo/Research/SimulationCode/FDPS-master/src" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


