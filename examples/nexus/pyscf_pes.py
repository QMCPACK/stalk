#!/usr/bin/env $python_exe

from numpy import savetxt
$pyscfimport

$system

$calculation

savetxt('value.out', [[e_scf, 0.0]])