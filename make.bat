@echo off

call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022

rem ----------------------------------------------------------------------
rem ifort compile
rem ----------------------------------------------------------------------
ifort .\src\iric.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\common.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\output.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\ss_nu.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\cell2node.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\cip3d.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\diffusion.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\newgrd.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\non_advection.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\o3upwind.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\omgpsi.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\upwind3d.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\w12cal.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\bound.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort .\src\nays3dv.f90  /Qopenmp /nostandard-realloc-lhs /MD /c
ifort *.obj .\lib\iriclib.lib -o ".\install\nays3dv.exe"

del *.obj
del *.mod
