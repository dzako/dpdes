#!/usr/bin/perl
#script for synchronizing data with the "dist" package
use File::Copy;

print "content-type: text/html \n\n"; #The header
$path="../capd_dynsys_dPDEs/dPDEs/";
copy("SHTest.cpp", $path."SHTest.cpp") or die "File cannot be copied.";
copy("config.h", $path."config.h") or die "File cannot be copied.";
copy("Equations.h", $path."Equations.h") or die "File cannot be copied.";
copy("DiffInclTest.cpp", $path."DiffInclTest.cpp") or die "File cannot be copied.";
copy("FFTTests.cpp", $path."FFTTests.cpp") or die "File cannot be copied.";
copy("PolyBdTests.cpp", $path."PolyBdTests.cpp") or die "File cannot be copied.";
copy("IndicesTests.cpp", $path."IndicesTests.cpp") or die "File cannot be copied.";
copy("KSPOTest.cpp", $path."KSPOTest.cpp") or die "File cannot be copied.";
copy("KSUnstPOTest.cpp", $path."KSUnstPOTest.cpp") or die "File cannot be copied.";
copy("BurgersFPTest.cpp", $path."BurgersFPTest.cpp") or die "File cannot be copied.";
copy("DiffInclusion.h", $path."DiffInclusion.h") or die "File cannot be copied.";
copy("DiffInclusion.hpp", $path."DiffInclusion.hpp") or die "File cannot be copied.";
copy("DPDEInclusionCW.h", $path."DPDEInclusionCW.h") or die "File cannot be copied.";
copy("PolyBd.h", $path."PolyBd.h") or die "File cannot be copied.";
copy("ComplexPolyBdJetOptimized.h", $path."ComplexPolyBdJetOptimized.h") or die "File cannot be copied.";
copy("FFT.h", $path."FFT.h") or die "File cannot be copied.";
copy("Index.h", $path."Index.h") or die "File cannot be copied.";
copy("DFTGrid.h", $path."DFTGrid.h") or die "File cannot be copied.";
copy("Real.h", $path."Real.h") or die "File cannot be copied.";
copy("DPDEContainer.h", $path."DPDEContainer.h") or die "File cannot be copied.";
copy("FirstOrderJet.h", $path."FirstOrderJet.h") or die "File cannot be copied.";
copy("FFTDynSys.h", $path."FFTDynSys.h") or die "File cannot be copied.";
copy("Even.h", $path."Even.h") or die "File cannot be copied.";
copy("Odd.h", $path."Odd.h") or die "File cannot be copied.";
copy("InclRect2Set.h", $path."InclRect2Set.h") or die "File cannot be copied.";
copy("ComplexScalar.h", $path."ComplexScalar.h") or die "File cannot be copied.";
copy("dissipative_enclosure.h", $path."dissipative_enclosure.h") or die "File cannot be copied.";
copy("tail2.h", $path."tail2.h") or die "File cannot be copied.";
copy("SetBlowUpException.h", $path."SetBlowUpException.h") or die "File cannot be copied.";
copy("Pair.h", $path."Pair.h") or die "File cannot be copied.";
copy("norms.h", $path."norms.h") or die "File cannot be copied.";
copy("Coefficients.h", $path."Coefficients.h") or die "File cannot be copied.";
copy("InclRect2Set.h", $path."InclRect2Set.h") or die "File cannot be copied.";
copy("InclRect2Set.hpp", $path."InclRect2Set.hpp") or die "File cannot be copied.";
copy("README", $path."README") or die "File cannot be copied.";
