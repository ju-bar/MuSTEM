@echo off
set PGILIB=C:\Program Files\PGI\win64\19.4
set PGICUDA=C:\Program Files\PGI\win64\2019\cuda\10.1
set INCL="-I%PGILIB%\include" "-I%PGICUDALIB%\include"
set LIBDIRS="-Wl,/libpath:%PGILIB%\lib" "-Wl,/libpath:%PGICUDALIB%\lib\x64"
set LIBS=cufft.lib
set OPTS=-Mpreprocess -DGPU -Dsingle_precision -Bstatic -Mbackslash -Mfree -mp -Mcuda=cuda10.1,ccall -fast -tp=px-64 -Minform=warn -O3
set OBJDIR=x64\objgpu
set OUTDIR=x64\Release
set OUTNAME=muSTEM_GPU
set COMPILEDFILES=0
set MAKEALL=0
if "%1"=="all" set MAKEALL=1
if not exist %OBJDIR% (
	echo - creating folder %OBJDIR%
	mkdir %OBJDIR%
)
if not exist %OUTDIR% (
	echo - creating folder %OUTDIR%
	mkdir %OUTDIR%
)
echo Compiling and linking %OUTNAME%
REM COMPILATION OF ALL REQUIRED F90 FILES IN OBJDIR
for %%f in (
		m_precision.f90
		m_string.f90
		m_user_input.f90
		mod_global_variables.f90
		GPU_routines\mod_cufft.f90
		GPU_routines\mod_cuda_setup.f90
		GPU_routines\mod_cuda_array_library.f90
		mod_CUFFT_wrapper.f90
		mod_output.f90
		m_numerical_tools.f90
		m_crystallography.f90
		m_plasmon.f90
		m_multislice.f90
		m_lens.f90
		m_tilt.f90
		quadpack.f90
		m_electron.f90
		m_absorption.f90
		GPU_routines\mod_cuda_potential.f90
		m_potential.f90
		mod_Hn0.f90
		GPU_routines\mod_cuda_ms.f90
		MS_utilities.f90
		s_qep_tem.f90
		s_qep_stem.f90
		s_absorptive_tem.f90
		s_absorptive_stem.f90
		muSTEM.f90
		) do (
	set COMPILED=0
	call :pgicompile %%f COMPILED
	set /a COMPILEDFILES=COMPILEDFILES+COMPILED
)
echo - total files compiled: %COMPILEDFILES%
if exist %OUTDIR%\%OUTNAME%.exe (
	if %COMPILEDFILES% equ 0 (
		echo - no source changes, linking skipped
		echo Done
		exit /b
	)
)
REM LINKING ALL OBJECT FILES REQUIRED FOR MUSTEM
echo - linking %OUTDIR%\%OUTNAME%.exe
pgfortran %OPTS% %LIBDIRS% -o %OUTDIR%\%OUTNAME% %OBJDIR%\s_absorptive_stem.obj %OBJDIR%\s_absorptive_tem.obj %OBJDIR%\s_qep_stem.obj %OBJDIR%\s_qep_tem.obj %OBJDIR%\MS_utilities.obj %OBJDIR%\mod_cuda_ms.obj %OBJDIR%\m_plasmon.obj %OBJDIR%\mod_Hn0.obj %OBJDIR%\m_potential.obj %OBJDIR%\mod_cuda_potential.obj %OBJDIR%\m_absorption.obj %OBJDIR%\m_electron.obj %OBJDIR%\quadpack.obj %OBJDIR%\m_tilt.obj %OBJDIR%\m_lens.obj %OBJDIR%\m_multislice.obj %OBJDIR%\m_crystallography.obj %OBJDIR%\m_numerical_tools.obj %OBJDIR%\mod_output.obj %OBJDIR%\mod_CUFFT_wrapper.obj %OBJDIR%\mod_cuda_array_library.obj %OBJDIR%\mod_cuda_setup.obj %OBJDIR%\mod_cufft.obj %OBJDIR%\mod_global_variables.obj %OBJDIR%\m_user_input.obj %OBJDIR%\m_string.obj %OBJDIR%\m_precision.obj %OBJDIR%\muSTEM.obj %LIBS%
echo Done
@pause
exit /b

REM SUBROUTINE pgicompile -------------------------
REM - compiles the file by default
REM - skips compile if the object file is newer than the source file
:pgicompile
set SRCFILE=%1
set OBJFILE=%OBJDIR%\%~n1.obj
set DOCOMPILE=%MAKEALL%
if exist %OBJFILE% (
	xcopy /L /D /Y %SRCFILE% %OBJFILE%|findstr /B /C:"1 " && set DOCOMPILE=1
) else (
	set DOCOMPILE=1
)
if %DOCOMPILE% equ 1 (
	echo - compiling %SRCFILE%
	pgfortran %OPTS% %INCL% -module %OBJDIR% -o %OBJFILE% -c %SRCFILE%
	set %~2=1
)
exit /b