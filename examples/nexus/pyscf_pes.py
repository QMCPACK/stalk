#!/usr/bin/env $python_exe

from numpy import savetxt
$pyscfimport

$system

$calculation

savetxt('value.dat', [[e_scf, 0.0]])