#Corre el MO mcola
cp esmeco_real.dat esmeco.dat
rm MOM.ecm
rm MOM.mc2
rm MOM.mcm
rm MOM.eva
rm MOM.psv
./MOM -mcmc 100000 -mcsave 100
