CXX = g++
CXXFLAGS = -std=c++17 -O2
TARGET = main

$(TARGET): main.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) main.cpp

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
