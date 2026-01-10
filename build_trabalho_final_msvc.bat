@echo off
del a.exe
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvarsall.bat" x86_amd64
cl trabalho_final.cpp /Fe:a.exe /EHsc /std:c++17 user32.lib gdi32.lib
del trabalho_final.obj
