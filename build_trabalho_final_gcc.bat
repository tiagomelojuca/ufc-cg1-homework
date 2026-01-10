@echo off
del a.exe
g++ -Wno-deprecated-declarations -municode -mwindows .\trabalho_final.cpp -lgdi32
