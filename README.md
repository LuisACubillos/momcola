Modelo Operativo de merluza de cola con cambios decadales en el reclutamiento
===========

Este repositorio contiene el código del modelo operativo de merluza de cola frente a Chile. Curin (2014) utilizó el modelo operativo para evaluar el desempeño del modelo de evaluación de stock (MES) considerando escenarios futuros de cambio decadale en el reclutamiento de merluza de cola.

El código del MO es una modificación del código original del MES. El MES fue escrito y documentado por Payá y Canales (2012) y sometido a revisión en el [Taller de Trabajo de Evaluación de Stock de Merluza de Cola 2011 ](https://sites.google.com/site/chsaw2011/home) realizado entre el 4 y 8 de julio de 2011, Viña del Mar, Chile.


This repository contains the code of the operative model (OM) of Patagonian grenadier off Chile. Curin (2014) used the operative model to evaluate the performance of the stock assessment model (SAM) under future scenarios of decadal change in the recruitment of Patagonian grenadier.

The code of the OM is  a modification of the original code of the SAM. The SAM was written and documented by Payá and Canales (2012) and submitted for review in the [Chilean Hoki Stock Assessment Workshop 2011](https://sites.google.com/site/chsaw2011/home) carried out in Viña del Mar (Chile), from 4 to 11 July 2011.

# Aspectos técnicos

El programa es escrito en código C++ como template de [AD Model Builder](http://www.admb-project.org) (Fournier et al., 2012). La última versión del MO se puede obtener del repositorio en [github](https://github.com/LuisACubillos/momcola).

The software is writting in C++ code as template of [AD Model Builder](http://www.admb-project.org) (Fournier et al., 2012). The latest version of the OM may obtaining from the repository in [github](https://github.com/LuisACubillos/momcola). 

## Requerimientos

* Un compilador C++
* AD Model Builder v. 11 o superior

## Clonación

El código fue escrito para Mac OS, de tal manera se puede acceder desde la terminal:

	cd ~
	git clone https://github.com/luisacubillos/momcola
	cd momcola                                                                                                                             

En windows se debe modificar los llamados al system y los archivos batch que correspondan.

# Instrucciones

La carpeta omsrc contiene el código y archivos del modelo operativo:

* MOM.tpl
* MOM.dat
* MOM_op.ctl
* ProyectaFs.ctl, fija el numero de tasas de explotación a evaluar como multiplicadores de 0.18.
* esmeco_real.dat, datos del estimador. 

La carpeta samsrc contiene el código y archivos del estimador:

* esmeco.tpl
* esmeco.dat
* esmeco_op.ctl
* ProyectaFs.ctl, este archivo debe ser idéntico al que usará MOM.


1) Compilar el template **MOM.tpl**

2) El archivo proyectaFs.ctl contiene los multiplicadores de una tasa de explotación objetivo de 0.18, y contiene tres valores.

3) El estimador está en la carpeta samsrc. Se debe compilar el template esmeco.tpl, y trasladar el ejecutable y los archivos esmeco.dat, esmeco_op.ctl, y esmeco_real.dat a la carpeta donde se compiló MOM y donde esté el ejecutable MOM.

4) Ejecutar mcmc, por ejemplo: ./mom -mcmc 10000 -mcsave 100

5) Para activar el modelo operativo MOM: ./mom -mceval

6) Esperar a que el modelo termine un ciclo de proyección de 40 años con tres estrategias de explotación.

7) Para modificar un escenario de proyección se debe editar MOM_op.ctl, específicamente las siguientes líneas 60 a la 67. Una vez modificada se puede correr un nuevo escenario, pero se debe asegurar de copiar los archivos mcmc resultantes de la corrida anterior.



