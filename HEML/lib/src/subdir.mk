################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CipherGD.cpp \
../src/GD.cpp \
../src/CipherSVM.cpp\
../src/TrainSVM.cpp\
../src/HEML.cpp \
../src/TestGD.cpp 

OBJS += \
./src/CipherGD.o \
./src/GD.o \
./src/CipherSVM.o\
./src/TrainSVM.o\
./src/HEML.o \
./src/TestGD.o 

CPP_DEPS += \
./src/CipherGD.d \
./src/GD.d \
./src/CipherSVM.d\
./src/TrainSVM.d\
./src/HEML.d \
./src/TestGD.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I../../../HEANN-master/heaan/src -I/usr/local/include -I/home/saerom/HEAAN-master/heaan/src -O3 -pthread -c -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


