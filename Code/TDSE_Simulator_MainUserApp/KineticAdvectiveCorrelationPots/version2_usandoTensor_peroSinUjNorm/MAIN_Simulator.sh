#!/bin/bash

option=1
until [[ $option -eq Out ]]; do
    clear
    echo ""
    echo "  Welcome to the time dependent Schroedinger Equation solver by XOA"
    echo ""
    echo "  --------------------------------------------------------------"
    echo ""
    echo "  We dispose of several different algorithms for simulating the time evolution of the Wave Function numerically, choose the desired one :"
    echo ""
    echo " 0) Not Implemented - CN 1D tdependentV"
	echo " 1) Cranck Nicolson 1D with Absorbing Boundary Conditions TIME INDEPENDENT POTENTIALS"
	echo " 1.5) Same as option 1 but with an option to add code patches outside and inside main()"
	echo " 2) Cranck Nicolson 1D with Absorbing Boundary Conditions TIME DEPENDENT POTENTIALS"
	echo " 3) Cranck Nicolson 2D TIME INDEPENDENT POTENTIALS"
	echo " 4) Cranck Nicolson 2D TIME DEPENDENT POTENTIALS"
	echo " 5) XO algorithm nD ZERO order Taylor approximation with ABCs and Time Dependent Potentials NO MULTIPROCESSING"
	echo " 6) XO algorithm nD ZERO order Taylor approximation with ABCs and Time Dependent Potentials YES MULTIPROCESSING"
	echo " 7) ANIMATED XO algorithm 2D ZERO order Taylor approximation with ABCs and Time Dependent Potentials - DONT Get Tensor Product"
	echo " 8) ANIMATED XO algorithm 2D ZERO order Taylor approximation with ABCs and Time Dependent Potentials - Get Tensor Product"
	echo " 9) Dissipative Imaginary Time Evolution 1D CN CONSTANT DISSIPATION in time"
	echo " 9.5) Dissipative Imaginary Time Evolution 1D CN TIME DEPENDANT DISSIPATION in time"
	echo " 10) Dissipative Imaginary Time Evolution 2D CN CONSTANT dissipation in time "
	echo " 11) Dissipative Imaginary Time Evolution XO algorithm"
	echo " Out) Close Simulator"
	echo ""
	echo " ----------------------------------"
	echo "     Introduce an option (1/2/3/4/5/6) "
	echo ""
	read option

    case $option in
    
    1)
    
    psiIni=""
    potential=""
    nx=1
    xmin=1.0
    xmax=1.0
    dt=1.0
    tIt=1
    numTrajs=0
    customTrajs=" "
    numTrajs=0
    mass=1.0
    hbar=1.0
    gif=L
    speed=1
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space x please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of x"
    read potential
    echo "Introduce the number of space divisions for the grid, xmin and xmax"
    read nx xmin xmax
    echo "Introduce the time iteration number "
    read tIt
    echo " Introduce the time step"
    read dt
     echo "Introduce the mass of the particle in the considered units:"
    read mass
    echo "Introduce the Normalised Planck Constant to be used (eg 1 for atomic units)"
    read hbar
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    echo "Introduce the number of time steps to jump between frames of the animation (if 1, every frame will be plotted, if 5, only one every 5 time steps will be plotted)"
    read speed
    
    echo "Do you want to plot Bohmian Trajectories? Options:"
    echo " N) no trajectories will be plotted"
    echo " R) randomly chosen trajectories will be plotted according to the initial |WF|**2"
    echo " C) you will choose the initial trajectory positions"
    read trajOption
    if [[ $trajOption == *"R"* || $trajOption == *"C"* ]]; then
        echo "Choose the number of trajectories:"
        read numTrajs
    fi
    if [[ $trajOption == *"C"* ]]; then
        echo "Introduce as one, space, another, space, another... the initial positions as numbers over the x axis between " $xmin " and " $xmax
        read customTrajs
    fi
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " $psiIni 
    echo " V(x) = " $potential 
    echo " nx = " $nx "; xmin = " $xmin "; xmax = " $xmax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo " hbar = " $habr "; mass = " $mass
    echo ""
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    codeOutsideMain=" "
    codeInsideMain= " "
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    ./EXE_codeFileGenerator_1D_CN "$psiIni" "$potential" $nx $xmin $xmax $dt $tIt $mass $hbar $option "$codeOutsideMain" "$codeInsideMain"
    
    g++ -Wall -O CODE_simulator_1D_CN_ABC_tINDEP.cpp -o EXE_simulator_1D_CN_ABC_tINDEP
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_1D_CN_ABC_tINDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    echo " Exporting plot data..."
    path=$(pwd)/DATA_rawSimulationData_1D_CN.txt
    ./EXE_prepareDataToPlot_1D_CN $path $trajOption $numTrajs "$customTrajs"
    EXPORT_TIME=$SECONDS
    echo " Done!" $(($EXPORT_TIME - $CALCULATION_TIME)) " seconds required!"
    
    echo ""
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."

    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  

        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tIndepPot_$currentDateTime nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; t=t+$speed;}"
        else
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tIndepPot_$currentDateTime nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' axes x1y1 w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' axes x1y1 w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories' axes x1y1 w p pt 6, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; t=t+$speed;}"
        fi
    
    
    else
    
        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; set terminal wxt size 1000,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; pause 0.1; t=t+$speed;}"
        else
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; set terminal wxt size 1000,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' axes x1y1 w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' axes x1y1 w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories' axes x1y1 w p pt 6, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; pause 0.1; t=t+$speed;}"
        fi
    fi
    
    
	;;
	
	1.5)
    
    psiIni=""
    potential=""
    nx=1
    xmin=1.0
    xmax=1.0
    dt=1.0
    tIt=1
    numTrajs=0
    customTrajs=" "
    numTrajs=0
    mass=1.0
    hbar=1.0
    gif=L
    speed=1
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space x please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of x"
    read potential
    echo "Introduce the number of space divisions for the grid, xmin and xmax"
    read nx xmin xmax
    echo "Introduce the time iteration number "
    read tIt
    echo " Introduce the time step"
    read dt
     echo "Introduce the mass of the particle in the considered units:"
    read mass
    echo "Introduce the Normalised Planck Constant to be used (eg 1 for atomic units)"
    read hbar
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    echo "Introduce the number of time steps to jump between frames of the animation (if 1, every frame will be plotted, if 5, only one every 5 time steps will be plotted)"
    read speed
    
    echo "Do you want to plot Bohmian Trajectories? Options:"
    echo " N) no trajectories will be plotted"
    echo " R) randomly chosen trajectories will be plotted according to the initial |WF|**2"
    echo " C) you will choose the initial trajectory positions"
    read trajOption
    if [[ $trajOption == *"R"* || $trajOption == *"C"* ]]; then
        echo "Choose the number of trajectories:"
        read numTrajs
    fi
    if [[ $trajOption == *"C"* ]]; then
        echo "Introduce as one, space, another, space, another... the initial positions as numbers over the x axis between " $xmin " and " $xmax
        read customTrajs
    fi
    
    echo " Introduce additional C++ code for the space outside the main function "
    read codeOutsideMain
    
    echo "Introduce additional C++ code to copy inside the main function after the space-time grid variable declaration"
    read codeInsideMain
    
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " $psiIni 
    echo " V(x) = " $potential 
    echo " nx = " $nx "; xmin = " $xmin "; xmax = " $xmax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo " hbar = " $habr "; mass = " $mass
    echo ""
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    ./EXE_codeFileGenerator_1D_CN "$psiIni" "$potential" $nx $xmin $xmax $dt $tIt $mass $hbar $option "$codeOutsideMain" "$codeInsideMain"
    
    g++ -Wall -O CODE_simulator_1D_CN_ABC_tINDEP.cpp -o EXE_simulator_1D_CN_ABC_tINDEP
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_1D_CN_ABC_tINDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    echo " Exporting plot data..."
    path=$(pwd)/DATA_rawSimulationData_1D_CN.txt
    ./EXE_prepareDataToPlot_1D_CN $path $trajOption $numTrajs "$customTrajs"
    EXPORT_TIME=$SECONDS
    echo " Done!" $(($EXPORT_TIME - $CALCULATION_TIME)) " seconds required!"
    
    echo ""
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."

    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  

        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tIndepPot_$currentDateTime nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; t=t+$speed;}"
        else
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tIndepPot_$currentDateTime nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' axes x1y1 w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' axes x1y1 w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories' axes x1y1 w p pt 6, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; t=t+$speed;}"
        fi
    
    
    else
    
        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; set terminal wxt size 1000,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; pause 0.1; t=t+$speed;}"
        else
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; set terminal wxt size 1000,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' axes x1y1 w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' axes x1y1 w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories' axes x1y1 w p pt 6, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; pause 0.1; t=t+$speed;}"
        fi
    fi
    
    
	;;

	2)
    
    psiIni=""
    potential=""
    nx=1
    xmin=1.0
    xmax=1.0
    dt=1.0
    tIt=1
    numTrajs=0
    mass=1.0
    hbar=1.0

    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space x please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of x and t"
    read potential
    echo "Introduce the number of space divisions for the grid, xmin and xmax"
    read nx xmin xmax
    echo "Introduce the time iteration number "
    read tIt
    echo " Introduce the time step"
    read dt
    echo "Introduce the mass of the particle in the considered units:"
    read mass
    echo "Introduce the Normalised Planck Constant to be used (eg 1 for atomic units)"
    read hbar
    echo ""
    gif=L
    speed=1
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    echo "Introduce the number of time steps to jump between frames of the animation (if 1, every frame will be plotted, if 5, only one every 5 time steps will be plotted)"
    read speed
    
     
    echo "Do you want to plot Bohmian Trajectories? Options:"
    echo " N) no trajectories will be plotted"
    echo " R) randomly chosen trajectories will be plotted according to the initial |WF|**2"
    echo " C) you will choose the initial trajectory positions"
    read trajOption
    if [[ $trajOption == *"R"* || $trajOption == *"C"* ]]; then
        echo "Choose the number of trajectories:"
        read numTrajs
    fi
    if [[ $trajOption == *"C"* ]]; then
        echo "Introduce as one, space, another, space, another... the initial positions as numbers over the x axis between " $xmin " and " $xmax
        read customTrajs
    fi
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " $psiIni 
    echo " V(x) = " $potential 
    echo " nx = " $nx "; xmin = " $xmin "; xmax = " $xmax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo " hbar = " $habr "; mass = " $mass
    echo ""
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    ./EXE_codeFileGenerator_1D_CN "$psiIni" "$potential" $nx $xmin $xmax $dt $tIt $mass $hbar $option
    
    g++ -Wall -O CODE_simulator_1D_CN_ABC_tDEP.cpp -o EXE_simulator_1D_CN_ABC_tDEP
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_1D_CN_ABC_tDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    echo " Exporting plot data..."
    path=$(pwd)/DATA_rawSimulationData_1D_CN.txt
    ./EXE_prepareDataToPlot_1D_CN $path $trajOption $numTrajs "$customTrajs" 
    EXPORT_TIME=$SECONDS
    echo " Done!" $(($EXPORT_TIME - $CALCULATION_TIME)) " seconds required!"
    echo ""
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."
    
    #TODO POTENTIAL TIME DEPENDANTA PLOTIEUUU!
    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`        
        echo " Generating Animation GIF..."  
        
        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tDepPot_$currentDateTime nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set yrange [-0.8:0.8]; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' w l; t=t+$speed;}"
        else
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tDepPot_$currentDateTime nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set yrange [-0.8:0.8]; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories'; t=t+$speed;}"
        fi
        
        
    else

        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set yrange [-0.8:0.8]; set xlabel 'Position x'; set terminal wxt size 900,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' w l; pause 0.1; t=t+$speed;}"
        else
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set yrange [-0.8:0.8]; set xlabel 'Position x'; set terminal wxt size 900,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories'; pause 0.1; t=t+$speed;}"
        fi
    fi
    
    ;;
	
	3)
	
	psiIni=""
    potential=""
    nx1=1
    nx2=1
    x1min=1.0
    x1max=1.0
    x2min=1.0
    x2max=1.0
    dt=1.0
    tIt=1
    plotOpt=1
    numTrajs=0
    mass1=1.0
    mass2=1.0
    hbar=1.0
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space variables x1 and x2 please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of x1 and x2"
    read potential
    echo "Introduce the number of space divisions for the grid in x1 and x2"
    read nx1 nx2
    echo "Introduce x1min x1max x2min x2max"
    read x1min x1max x2min x2max
    echo "Introduce the time iteration number "
    read tIt
    echo " Introduce the time step"
    read dt
    echo " Introduce the mass at x1 and mass at x2 as m1 m2"
    read mass1 mass2
    echo " Introduce the value of hbar"
    read hbar
    echo "Introduce the plot modality 1) 2D heat map or 2) 3D Probability surface"
    read plotOpt
    echo ""
    gif=L
    speed=1
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    echo "Introduce the number of time steps to jump between frames of the animation (if 1, every frame will be plotted, if 5, only one every 5 time steps will be plotted)"
    read speed
    
    echo "Do you want to plot Bohmian Trajectories? Options:"
    echo " N) no trajectories will be plotted"
    echo " R) randomly chosen trajectories will be plotted according to the initial |WF|**2"
    echo " C) you will choose the initial trajectory positions"
    read trajOption
    if [[ $trajOption == *"R"* || $trajOption == *"C"* ]]; then
        echo "Choose the number of trajectories:"
        read numTrajs
    fi
    if [[ $trajOption == *"C"* ]]; then
        echo "Introduce as a1, space, b1, space, a2, space, b2,... the initial positions as numbers over the x1 axis between " $x1min " and " $x1max  "and the x2 axis between "$x2min " and " $x2max
        read customTrajs
    fi
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " "$psiIni"
    echo " V(x) = " "$potential" 
    echo " nx1 = " $nx1 "; nx2 = " $nx2 
    echo " x1min = " $x1min "; x1max = " $x1max "; x2min = " $x2min "; x2max = " $x2max
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo ""
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    ./EXE_codeFileGenerator_2D_CN "$psiIni" "$potential" $nx1 $nx2 $x1min $x1max $x2min $x2max $dt $tIt $mass1 $mass2 $hbar $option
    g++ -Wall -O CODE_simulator_2D_CN_tINDEP.cpp -o EXE_simulator_2D_CN_tINDEP
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_2D_CN_tINDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    echo " Exporting plot data..."
    path=$(pwd)/DATA_rawSimulationData_2D_CN.txt
    ./EXE_prepareDataToPlot_2D_CN $path $trajOption $numTrajs "$customTrajs" 
    EXPORT_TIME=$SECONDS
    echo " Done!" $(($EXPORT_TIME - $CALCULATION_TIME)) " seconds required!"
    echo ""
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."
    
    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  

        if [[ $trajOption == *"N"* ]]; then
        
        if [[ $plotOpt == 1 ]]; then
            gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN2D_tIndepPot_$currentDateTime nx1$nx1 nx2$nx2 x1min$x1min x2min$x2min x1max$x1max x2max$x2max dt$dt tIt$tIt mass1$mass1 mass2$mass2 hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x1'; set ylabel 'Position x2'; set contour base; set cntrparam levels 3; set palette rgbformulae 33,13,10; set pm3d map; set key at 25,15 title 'Wave Function Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1); t=t+$speed;}"
        elif [ $plotOpt == 2 ]; then
            gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN2D_tIndepPot_$currentDateTime nx1$nx1 nx2$nx2 x1min$x1min x2min$x2min x1max$x1max x2max$x2max dt$dt tIt$tIt mass1$mass1 mass2$mass2 hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set zrange [0:0.13]; set xlabel 'Position x1'; set ylabel 'Position x2'; unset hidden3d; set pm3d at bs; set palette rgbformulae 33,13,10; set contour; unset surface; set key title 'WaveFunction Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1) w l; t=t+$speed;}"
        fi

        else
        
            gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN2D_tIndepPot_$currentDateTime nx1$nx1 nx2$nx2 x1min$x1min x2min$x2min x1max$x1max x2max$x2max dt$dt tIt$tIt mass1$mass1 mass2$mass2 hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x1'; set ylabel 'Position x2'; set palette rgbformulae 33,13,10; set pm3d map; set key at -4,10 title 'Wave Function Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1), 'DATA_trajectoriesToPlot_2D_CN.txt' index t title 'Trajectories' w points lt 5 pt 7 ps 1 lc rgb 'black'; t=t+$speed;}"
        
        
        fi
        
    else
    
        if [[ $trajOption == *"N"* ]]; then
        
        if [[ $plotOpt == 1 ]]; then
            gnuplot -e "t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x1'; set ylabel 'Position x2'; set contour base; set cntrparam levels 3; set palette rgbformulae 33,13,10; set pm3d map; set terminal wxt size 900,900; set key at 25,15 title 'Wave Function Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1); pause 0.00001; t=t+$speed;}"
        elif [ $plotOpt == 2 ]; then
            gnuplot -e "t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set zrange [0:0.13]; set xlabel 'Position x1'; set ylabel 'Position x2'; unset hidden3d; set pm3d at bs; set palette rgbformulae 33,13,10; set contour; unset surface;  set key title 'WaveFunction Probability Density'; set terminal wxt size 900,900; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1) w l; pause 0.00001; t=t+$speed;}"
        fi

        else
        
            gnuplot -e "t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x1'; set ylabel 'Position x2'; set palette rgbformulae 33,13,10; set pm3d map; set terminal wxt size 900,900; set key at -4,10 title 'Wave Function Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1), 'DATA_trajectoriesToPlot_2D_CN.txt' index t title 'Trajectories' w points lt 5 pt 7 ps 1 lc rgb 'black' ; pause 0.00001; t=t+$speed;}"
        
        
        fi
    fi
    
	;;
	
	5)
	psiIni=""
    potential=""
    qDivs=""
    qmin=""
    qmax=""
    dt=1.0
    tIt=1
    plotOpt=1
    numTrajs=0
    dimNum=1
    mass=""
    hbar=1.0
    potentialPlotFineness=1.0
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce n, the number of spatial dimensions of the problem"
    read dimNum
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space variables q[0], q[1],...,q[n] please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of q[0],..,q[n]"
    read potential
    echo "Introduce the mass for each dimension separated by commas m[0], m[1], ... m[n]"
    read mass
    echo "Introduce the number of space divisions for the grid in each dimension  separated by commas q[0], q[1], ... q[n]"
    read qDivs
    echo "Introduce the lower grid bounds qmin[0], qmin[1], ... qmin[n]"
    read qmin
    echo "Introduce the upper grid bounds qmax[0], qmax[1], ... qmax[n]"
    read qmax
    echo "Introduce the time iteration number "
    read tIt
    echo "Introduce the time step"
    read dt
    echo "Introduce the number of trajectories to compute"
    read numTrajs
    echo "Introduce hbar"
    read hbar
    echo "If the system is 2D you can plot the potential profile as a heat map. Introduce the plot modality 1) dont plot potential or 2) plot potential field"
    read plotOpt
    echo ""
    gif=L
    speed=1
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " $psiIni 
    echo " V(x) = " $potential 
    echo " qDivs = " $qDivs  
    echo " qmin = " $qmin
    echo " qmax = " $qmax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo ""
    
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    
    ./EXE_codeFileGenerator_nD_XO_ZERO_CN_ABC_tDEP $dimNum "$psiIni" "$potential" "$mass" "$qDivs" "$qmin" "$qmax" $dt $tIt $option $numTrajs $potentialPlotFineness $hbar

    g++ -Wall -O CODE_simulator_nD_XO_ZERO_CN_ABC_tDEP.cpp -o EXE_simulator_nD_XO_ZERO_CN_ABC_tDEP -lm
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_nD_XO_ZERO_CN_ABC_tDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    #echo " Exporting plot data..."
    #path=$(pwd)/DATA_rawSimulationData_2D_CN.txt
    #./EXE_prepareDataToPlot_2D_CN $path $trajOption $numTrajs "$customTrajs" 
    #EXPORT_TIME=$SECONDS
    #echo " Done!" $(($EXPORT_TIME - $CALCULATION_TIME)) " seconds required!"
    #echo ""
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."
    
    x1min=$(echo $qmin| cut -d',' -f 1)
    x2min=$(echo $qmin| cut -d',' -f 2)
    x1max=$(echo $qmax| cut -d',' -f 1)
    x2max=$(echo $qmax| cut -d',' -f 2)
    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel XO_$currentDateTime dim$dimNum qDivs($qDivs) qmin($qmin) qmax($qmax) numTrajs($numTrajs) dt$dt tIt$tIt mass($mass) hbar$hbar.gif'; t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; set pm3d map; splot 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::1::1 title 'Initial Coordinates' w points lt 5 pt 6 ps 2 lc rgb 'red', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::($tIt-1)::($tIt-1) title 'Final Coordinates' w points lt 5 pt 6 ps 2 lc rgb 'blue';"
    
    
    else
        gnuplot --persist -e "t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; set terminal wxt size 900,900; set pm3d map; splot 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::1::1 title 'Initial Coordinates' w points lt 5 pt 6 ps 2 lc rgb 'red', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::($tIt-1)::($tIt-1) title 'Final Coordinates' w points lt 5 pt 6 ps 2 lc rgb 'blue'; pause 20; reread;"
    
    fi
	
	;;
	
	7)
	psiIni=""
    potential=""
    qDivs=""
    qmin=""
    qmax=""
    dt=1.0
    tIt=1
    plotOpt=1
    numTrajs=0
    dimNum=1
    mass=""
    speedUp=1
    potentialPlotFineness=0.02
    hbar=1.0
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce n, the number of spatial dimensions of the problem"
    read dimNum
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space variables q[0], q[1],...,q[n] please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of q[0],..,q[n]"
    read potential
    echo "Introduce the mass for each dimension separated by commas m[0], m[1], ... m[n]"
    read mass
    echo "Introduce the number of space divisions for the grid in each dimension  separated by commas q[0], q[1], ... q[n]"
    read qDivs
    echo "Introduce the lower grid bounds qmin[0], qmin[1], ... qmin[n]"
    read qmin
    echo "Introduce the upper grid bounds qmax[0], qmax[1], ... qmax[n]"
    read qmax
    echo "Introduce the time iteration number "
    read tIt
    echo "Introduce the time step"
    read dt
    echo "Introduce the number of trajectories to compute"
    read numTrajs
    echo "Introduce hbar"
    read hbar
    echo "Introduce the fineness with which you wish to plot the potential energy field  (0, 0.5). The lower the better, a good deal time-quality might be around 0.03"
    read potentialPlotFineness
    echo "The animation can be accelerated if you wish. Enter a speedUp factor (1/2/3/4/5/...10)"
    read speedUp
    echo ""
    gif=L
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " $psiIni 
    echo " V(x) = " $potential 
    echo " qDivs = " $qDivs  
    echo " qmin = " $qmin
    echo " qmax = " $qmax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo ""
    
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    
    ./EXE_codeFileGenerator_nD_XO_ZERO_CN_ABC_tDEP $dimNum "$psiIni" "$potential" "$mass" "$qDivs" "$qmin" "$qmax" $dt $tIt $option $numTrajs $potentialPlotFineness $hbar
    
    g++ -Wall -O CODE_simulator_nD_XO_ZERO_CN_ABC_tDEP.cpp -o EXE_simulator_nD_XO_ZERO_CN_ABC_tDEP -lm
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_nD_XO_ZERO_CN_ABC_tDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."
    
    x1min=$(echo $qmin| cut -d',' -f 1)
    x2min=$(echo $qmin| cut -d',' -f 2)
    x1max=$(echo $qmax| cut -d',' -f 1)
    x2max=$(echo $qmax| cut -d',' -f 2)
    divsx1=$(echo $qDivs| cut -d',' -f 1)
    divsx2=$(echo $qDivs| cut -d',' -f 2)
    
    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  
        
        
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel XO_$currentDateTime dim$dimNum qDivs($qDivs) qmin($qmin) qmax($qmax) numTrajs($numTrajs) dt$dt tIt$tIt mass($mass) hbar$hbar potentialPlotFineness$potentialPlotFineness speed$speedUp.gif'; t=0; countTrajs=0; tmax=2*$tIt*$numTrajs; dx1=($x1max-($x1min))/$divsx1.0; dx2=($x2max-($x2min))/$divsx2.0; set multiplot; array posx1[$divsx1+1]; array posx2[$divsx2+1]; do for [i=1:($divsx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($divsx2+1)] { posx2[i] = $x2min + i*dx2 }; while(t<tmax){ 
        set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.025]; set xlabel 'Position q1'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index t using (posx1[\$0+1]):1 title 'q1 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 1:3 title 'q1 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$x2min:$x2max]; set yrange [0:0.025]; set xlabel 'Position q2'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t+1) using (posx2[\$0+1]):1 title 'q2 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 2:3 title 'q2 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; if ((t-(2*$tIt*countTrajs)) >= 2*$tIt){ countTrajs=countTrajs+1; }; set palette rgbformulae 33,13,10; set pm3d map; splot 'DATA_potentialToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index ((t-(2*$tIt*countTrajs))/2) title 'Potential Energy Map', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::1::(t/2+1) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red'; 
        t=t+2*$speedUp; } unset multiplot;"
        
    else
        
        gnuplot -e "t=0; countTrajs=0; tmax=2*$tIt*$numTrajs; dx1=($x1max-($x1min))/$divsx1.0; dx2=($x2max-($x2min))/$divsx2.0;  set terminal wxt size 1800,900; set multiplot; array posx1[$divsx1+1]; array posx2[$divsx2+1]; do for [i=1:($divsx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($divsx2+1)] { posx2[i] = $x2min + i*dx2 }; while(t<tmax){ 
        set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.025]; set xlabel 'Position q1'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index t using (posx1[\$0+1]):1 title 'q1 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 1:3 title 'q1 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$x2min:$x2max]; set yrange [0:0.025]; set xlabel 'Position q2'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t+1) using (posx2[\$0+1]):1 title 'q2 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 2:3 title 'q2 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; if ((t-(2*$tIt*countTrajs)) >= 2*$tIt){ countTrajs=countTrajs+1; }; set palette rgbformulae 33,13,10; set pm3d map; splot 'DATA_potentialToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index ((t-(2*$tIt*countTrajs))/2) title 'Potential Energy Map', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::1::(t/2+1) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red'; 
        t=t+2*$speedUp; } unset multiplot; pause 3;"
	fi
	;;
	
	8)
	psiIni=""
    potential=""
    qDivs=""
    qmin=""
    qmax=""
    dt=1.0
    tIt=1
    plotOpt=1
    numTrajs=0
    dimNum=1
    mass=""
    speedUp=1
    potentialPlotFineness=0.02
    hbar=1.0
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce n, the number of spatial dimensions of the problem"
    read dimNum
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space variables q[0], q[1],...,q[n] please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of q[0],..,q[n]"
    read potential
    echo "Introduce the mass for each dimension separated by commas m[0], m[1], ... m[n]"
    read mass
    echo "Introduce the number of space divisions for the grid in each dimension  separated by commas q[0], q[1], ... q[n]"
    read qDivs
    echo "Introduce the lower grid bounds qmin[0], qmin[1], ... qmin[n]"
    read qmin
    echo "Introduce the upper grid bounds qmax[0], qmax[1], ... qmax[n]"
    read qmax
    echo "Introduce the time iteration number "
    read tIt
    echo "Introduce the time step"
    read dt
    echo "Introduce the number of trajectories to compute"
    read numTrajs
    echo "Introduce hbar"
    read hbar
    echo "Introduce the fineness with which you wish to plot the potential energy field and the tensor product (0, 0.5). The lower the better, a good deal time-quality might be around 0.03"
    read potentialPlotFineness
    echo "The animation can be accelerated if you wish. Enter a speedUp factor (1/2/3/4/5/...10)"
    read speedUp
    echo ""
    gif=L
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " $psiIni 
    echo " V(x) = " $potential 
    echo " qDivs = " $qDivs  
    echo " qmin = " $qmin
    echo " qmax = " $qmax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo ""
    
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    
    ./EXE_codeFileGenerator_nD_XO_ZERO_CN_ABC_tDEP $dimNum "$psiIni" "$potential" "$mass" "$qDivs" "$qmin" "$qmax" $dt $tIt $option $numTrajs $potentialPlotFineness $hbar

    g++ -Wall -O CODE_simulator_nD_XO_ZERO_CN_ABC_tDEP.cpp -o EXE_simulator_nD_XO_ZERO_CN_ABC_tDEP -lm
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_nD_XO_ZERO_CN_ABC_tDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."
    
    x1min=$(echo $qmin| cut -d',' -f 1)
    x2min=$(echo $qmin| cut -d',' -f 2)
    x1max=$(echo $qmax| cut -d',' -f 1)
    x2max=$(echo $qmax| cut -d',' -f 2)
    divsx1=$(echo $qDivs| cut -d',' -f 1)
    divsx2=$(echo $qDivs| cut -d',' -f 2)
    
    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  
        
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel XO_TensorProduct_$currentDateTime dim$dimNum qDivs($qDivs) qmin($qmin) qmax($qmax) numTrajs($numTrajs) dt$dt tIt$tIt mass($mass) hbar$hbar potentialPlotFineness$potentialPlotFineness speed$speedUp.gif'; t=0; countTrajs=0; tmax=2*$tIt*$numTrajs; dx1=($x1max-($x1min))/$divsx1.0; dx2=($x2max-($x2min))/$divsx2.0; array posx1[$divsx1+1]; array posx2[$divsx2+1]; do for [i=1:($divsx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($divsx2+1)] { posx2[i] = $x2min + i*dx2 }; while(t<tmax){  set multiplot; set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.025]; set xlabel 'Position q1'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index t using (posx1[\$0+1]):1 title 'q1 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 1:3 title 'q1 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$x2min:$x2max]; set yrange [0:0.025]; set xlabel 'Position q2'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t+1) using (posx2[\$0+1]):1 title 'q2 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 2:3 title 'q2 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; if ((t-(2*$tIt*countTrajs)) >= 2*$tIt){ countTrajs=countTrajs+1; }; set pm3d map; set palette rgbformulae 33,13,10; splot 'DATA_potentialToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index ((t-(2*$tIt*countTrajs))/2) title 'Potential Energy Map', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::1::(t/2+1) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';
        set origin 0.49, 0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; set palette rgbformulae 33,13,10; splot 'DATA_tensorProduct_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t/2);
        t=t+2*$speedUp; unset multiplot; }"
        echo "Done!"
    else
        gnuplot -e "t=0; countTrajs=0; tmax=2*$tIt*$numTrajs; dx1=($x1max-($x1min))/$divsx1.0; dx2=($x2max-($x2min))/$divsx2.0;  set terminal wxt size 1850,970; set multiplot; set pm3d map; array posx1[$divsx1+1]; array posx2[$divsx2+1]; do for [i=1:($divsx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($divsx2+1)] { posx2[i] = $x2min + i*dx2 }; while(t<tmax){ 
        set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.025]; set xlabel 'Position q1'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index t using (posx1[\$0+1]):1 title 'q1 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 1:3 title 'q1 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$x2min:$x2max]; set yrange [0:0.025]; set xlabel 'Position q2'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t+1) using (posx2[\$0+1]):1 title 'q2 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 2:3 title 'q2 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; if ((t-(2*$tIt*countTrajs)) >= 2*$tIt){ countTrajs=countTrajs+1; }; set palette rgbformulae 33,13,10; splot 'DATA_potentialToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index ((t-(2*$tIt*countTrajs))/2) title 'Potential Energy Map', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::1::(t/2+1) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';
        set origin 0.49, 0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; set palette rgbformulae 33,13,10; splot 'DATA_tensorProduct_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t/2);
        t=t+2*$speedUp; } unset multiplot; pause 3;"
            
    	#unset view; set pm3d;
    fi
	;;
	
	9)
    
    psiIni=""
    potential=""
    nx=1
    xmin=1.0
    xmax=1.0
    dt=1.0
    tIt=1
    numTrajs=0
    customTrajs=" "
    numTrajs=0
    mass=1.0
    hbar=1.0
    gif=L
    speed=1
    substituteFor_i=1.0
    codeOutsideMain=""
    codeInsideMain= ""
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space x please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of x"
    read potential
    echo "Introduce the number of space divisions for the grid, xmin and xmax"
    read nx xmin xmax
    echo "Introduce the time iteration number "
    read tIt
    echo " Introduce the time step"
    read dt
     echo "Introduce the mass of the particle in the considered units:"
    read mass
    echo "Introduce the Normalised Planck Constant to be used (eg 1 for atomic units)"
    read hbar
    echo "Introduce the desired substitution for the present i in the Schrodinger Equation for the dissipative time evolution (for ex. use 1 for a full imaginary time evolution, use J-(J -1)*eps with eps in (0,1) for a gradual linear interpolation between dissipative and standard time evolution...)"
    read substituteFor_i
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    echo "Introduce the number of time steps to jump between frames of the animation (if 1, every frame will be plotted, if 5, only one every 5 time steps will be plotted)"
    read speed
    
    echo "Do you want to plot Bohmian Trajectories? Options:"
    echo " N) no trajectories will be plotted"
    echo " R) randomly chosen trajectories will be plotted according to the initial |WF|**2"
    echo " C) you will choose the initial trajectory positions"
    read trajOption
    if [[ $trajOption == *"R"* || $trajOption == *"C"* ]]; then
        echo "Choose the number of trajectories:"
        read numTrajs
    fi
    if [[ $trajOption == *"C"* ]]; then
        echo "Introduce as one, space, another, space, another... the initial positions as numbers over the x axis between " $xmin " and " $xmax
        read customTrajs
    fi
    
    echo " Introduce additional C++ code for the space outside the main function "
    read codeOutsideMain
    
    echo "Introduce additional C++ code to copy inside the main function after the space-time grid variable declaration"
    read codeInsideMain
    
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " $psiIni 
    echo " V(x) = " $potential 
    echo " nx = " $nx "; xmin = " $xmin "; xmax = " $xmax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo " hbar = " $habr "; mass = " $mass
    echo ""
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    ./EXE_codeFileGenerator_1D_CN "$psiIni" "$potential" $nx $xmin $xmax $dt $tIt $mass $hbar $option "$codeOutsideMain" "$codeInsideMain" "$substituteFor_i"
    
    g++ -Wall -O CODE_simulator_1D_CN_tINDEP_ImTimeEv.cpp -o EXE_simulator_1D_CN_tINDEP_ImTimeEv
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_1D_CN_tINDEP_ImTimeEv
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    echo " Exporting plot data..."
    path=$(pwd)/DATA_rawSimulationData_1D_CN.txt
    ./EXE_prepareDataToPlot_1D_CN $path $trajOption $numTrajs "$customTrajs"
    EXPORT_TIME=$SECONDS
    echo " Done!" $(($EXPORT_TIME - $CALCULATION_TIME)) " seconds required!"
    
    echo ""
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."

    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  

        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tIndepPot_ImTimeEvol_$currentDateTime substituteFor_i $substituteFor_i nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; t=t+$speed;}"
        else
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tIndepPot_ImTimeEvol_$currentDateTime substituteFor_i $substituteFor_i nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' axes x1y1 w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' axes x1y1 w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories' axes x1y1 w p pt 6, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; t=t+$speed;}"
        fi
    
    
    else
    
        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; set terminal wxt size 1000,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; pause 0.1; t=t+$speed;}"
        else
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set y2tics nomirror; set y2label 'Energy'; set ylabel 'WaveFunction'; set xlabel 'Position x'; set terminal wxt size 1000,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) axes x1y1 w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' axes x1y1 w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' axes x1y1 w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories' axes x1y1 w p pt 6, 'DATA_potentialData_1D_CN.txt' using 1:2 title 'Potential Energy U(x)' axes x1y2 w l; pause 0.1; t=t+$speed;}"
        fi
    fi
    
    
	;;
	
	9.5)
    
    psiIni=""
    potential=""
    nx=1
    xmin=1.0
    xmax=1.0
    dt=1.0
    tIt=1
    numTrajs=0
    customTrajs=" "
    numTrajs=0
    mass=1.0
    hbar=1.0
    gif=L
    speed=1
    substituteFor_i=1.0
    codeOutsideMain=" "
    codeInsideMain= " "
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space x please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of x"
    read potential
    echo "Introduce the number of space divisions for the grid, xmin and xmax"
    read nx xmin xmax
    echo "Introduce the time iteration number "
    read tIt
    echo " Introduce the time step"
    read dt
     echo "Introduce the mass of the particle in the considered units:"
    read mass
    echo "Introduce the Normalised Planck Constant to be used (eg 1 for atomic units)"
    read hbar
    echo "Introduce the desired substitution for the present i in the Schrodinger Equation for the dissipative time evolution (for ex. use 1 for a full imaginary time evolution, use J-(J -1)*eps with eps in (0,1) for a gradual linear interpolation between dissipative and standard time evolution...)"
    read substituteFor_i
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    echo "Introduce the number of time steps to jump between frames of the animation (if 1, every frame will be plotted, if 5, only one every 5 time steps will be plotted)"
    read speed
    
    echo "Do you want to plot Bohmian Trajectories? Options:"
    echo " N) no trajectories will be plotted"
    echo " R) randomly chosen trajectories will be plotted according to the initial |WF|**2"
    echo " C) you will choose the initial trajectory positions"
    read trajOption
    if [[ $trajOption == *"R"* || $trajOption == *"C"* ]]; then
        echo "Choose the number of trajectories:"
        read numTrajs
    fi
    if [[ $trajOption == *"C"* ]]; then
        echo "Introduce as one, space, another, space, another... the initial positions as numbers over the x axis between " $xmin " and " $xmax
        read customTrajs
    fi
    
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " $psiIni 
    echo " V(x) = " $potential 
    echo " nx = " $nx "; xmin = " $xmin "; xmax = " $xmax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo " hbar = " $habr "; mass = " $mass
    echo ""
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    if [[ $substituteFor_i != *"return"* ]]; then
        potential="return ${substituteFor_i}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    ./EXE_codeFileGenerator_1D_CN "$psiIni" "$potential" $nx $xmin $xmax $dt $tIt $mass $hbar "$option" "$codeOutsideMain" "$codeInsideMain" "$substituteFor_i"
    
    g++ -Wall -O CODE_simulator_1D_CN_tDEP_ImTimeEv.cpp -o EXE_simulator_1D_CN_tDEP_ImTimeEv
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_1D_CN_tDEP_ImTimeEv
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    echo " Exporting plot data..."
    path=$(pwd)/DATA_rawSimulationData_1D_CN.txt
    ./EXE_prepareDataToPlot_1D_CN $path $trajOption $numTrajs "$customTrajs"
    EXPORT_TIME=$SECONDS
    echo " Done!" $(($EXPORT_TIME - $CALCULATION_TIME)) " seconds required!"
    
    echo ""
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."

    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  

        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tDepPot_$currentDateTime substituteFor_i $substituteFor_i nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' w l; t=t+$speed;}"
        else
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN1D_tDepPot_$currentDateTime substituteFor_i $substituteFor_i nx$nx xmin$xmin xmax$xmax dt$dt tIt$tIt mass$mass hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set xlabel 'Position x'; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories'; t=t+$speed;}"
        fi
    
    
    else
    
        
        if [[ $trajOption == *"N"* ]]; then
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set xlabel 'Position x'; set terminal wxt size 900,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' w l; pause 0.1; t=t+$speed;}"
        else
        gnuplot -e "t=0; tmax=$tIt; set xrange [$xmin:$xmax]; set xlabel 'Position x'; set terminal wxt size 900,900; while(t<tmax){ plot 'DATA_dataToPlot_1D_CN.txt' index t using 1:2 title columnheader(1) w l, 'DATA_dataToPlot_1D_CN.txt' index t using 1:3 title 'Real(Psi)' w l,'DATA_dataToPlot_1D_CN.txt' index t using 1:4 title 'Imaginary(Psi)' w l,  'DATA_trajectoriesToPlot_1D_CN.txt' index t title 'Trajectories'; pause 0.1; t=t+$speed;}"
        fi
    fi
    
    
	;;
	
	10)
	
	psiIni=""
    potential=""
    nx1=1
    nx2=1
    x1min=1.0
    x1max=1.0
    x2min=1.0
    x2max=1.0
    dt=1.0
    tIt=1
    plotOpt=1
    numTrajs=0
    mass1=1.0
    mass2=1.0
    hbar=1.0
    substituteFor_i = ""
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space variables x1 and x2 please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of x1 and x2"
    read potential
    echo "Introduce the number of space divisions for the grid in x1 and x2"
    read nx1 nx2
    echo "Introduce x1min x1max x2min x2max"
    read x1min x1max x2min x2max
    echo "Introduce the time iteration number "
    read tIt
    echo " Introduce the time step"
    read dt
    echo " Introduce the mass at x1 and mass at x2 as m1 m2"
    read mass1 mass2
    echo " Introduce the value of hbar"
    read hbar
    echo "Introduce the desired substitution for the present i in the Schrodinger Equation for the dissipative time evolution (for ex. use 1 for a full imaginary time evolution, use J-(J -1)*eps with eps in (0,1) for a gradual linear interpolation between dissipative and standard time evolution...)"
    read substituteFor_i
    echo "Introduce the plot modality 1) 2D heat map or 2) 3D Probability surface"
    read plotOpt
    echo ""
    gif=L
    speed=1
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    echo "Introduce the number of time steps to jump between frames of the animation (if 1, every frame will be plotted, if 5, only one every 5 time steps will be plotted)"
    read speed
    
    echo "Do you want to plot Bohmian Trajectories? Options:"
    echo " N) no trajectories will be plotted"
    echo " R) randomly chosen trajectories will be plotted according to the initial |WF|**2"
    echo " C) you will choose the initial trajectory positions"
    read trajOption
    if [[ $trajOption == *"R"* || $trajOption == *"C"* ]]; then
        echo "Choose the number of trajectories:"
        read numTrajs
    fi
    if [[ $trajOption == *"C"* ]]; then
        echo "Introduce as a1, space, b1, space, a2, space, b2,... the initial positions as numbers over the x1 axis between " $x1min " and " $x1max  "and the x2 axis between "$x2min " and " $x2max
        read customTrajs
    fi
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " "$psiIni"
    echo " V(x) = " "$potential" 
    echo " nx1 = " $nx1 "; nx2 = " $nx2 
    echo " x1min = " $x1min "; x1max = " $x1max "; x2min = " $x2min "; x2max = " $x2max
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo ""
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    ./EXE_codeFileGenerator_2D_CN "$psiIni" "$potential" $nx1 $nx2 $x1min $x1max $x2min $x2max $dt $tIt $mass1 $mass2 $hbar $option "$substituteFor_i"
    g++ -Wall -O CODE_simulator_2D_CN_tINDEP_ImTimeEv.cpp -o EXE_simulator_2D_CN_tINDEP_ImTimeEv
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_2D_CN_tINDEP_ImTimeEv
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    echo " Exporting plot data..."
    path=$(pwd)/DATA_rawSimulationData_2D_CN.txt
    ./EXE_prepareDataToPlot_2D_CN $path $trajOption $numTrajs "$customTrajs" 
    EXPORT_TIME=$SECONDS
    echo " Done!" $(($EXPORT_TIME - $CALCULATION_TIME)) " seconds required!"
    echo ""
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."
    
    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  

        if [[ $trajOption == *"N"* ]]; then
        
        if [[ $plotOpt == 1 ]]; then
            gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN2D_ImagTimeEv_tIndepPot_$currentDateTime substituteFor_i $substituteFor_i nx1$nx1 nx2$nx2 x1min$x1min x2min$x2min x1max$x1max x2max$x2max dt$dt tIt$tIt mass1$mass1 mass2$mass2 hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x1'; set ylabel 'Position x2'; set contour base; set cntrparam levels 3; set palette rgbformulae 33,13,10; set pm3d map; set key at 25,15 title 'Wave Function Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1); t=t+$speed;}"
        elif [ $plotOpt == 2 ]; then
            gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN2D_ImagTimeEv_tIndepPot_$currentDateTime substituteFor_i $substituteFor_i nx1$nx1 nx2$nx2 x1min$x1min x2min$x2min x1max$x1max x2max$x2max dt$dt tIt$tIt mass1$mass1 mass2$mass2 hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set zrange [0:0.13]; set xlabel 'Position x1'; set ylabel 'Position x2'; unset hidden3d; set pm3d at bs; set palette rgbformulae 33,13,10; set contour; unset surface; set key title 'WaveFunction Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1) w l; t=t+$speed;}"
        fi

        else
        
            gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel CN2D_ImagTimeEv_tIndepPot_$currentDateTime substituteFor_i $substituteFor_i nx1$nx1 nx2$nx2 x1min$x1min x2min$x2min x1max$x1max x2max$x2max dt$dt tIt$tIt mass1$mass1 mass2$mass2 hbar$hbar trajs$trajOption speed$speed.gif'; t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x1'; set ylabel 'Position x2'; set palette rgbformulae 33,13,10; set pm3d map; set key at -4,10 title 'Wave Function Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1), 'DATA_trajectoriesToPlot_2D_CN.txt' index t title 'Trajectories' w points lt 5 pt 7 ps 1 lc rgb 'black'; t=t+$speed;}"
        
        
        fi
        
    else
    
        if [[ $trajOption == *"N"* ]]; then
        
        if [[ $plotOpt == 1 ]]; then
            gnuplot -e "t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x1'; set ylabel 'Position x2'; set contour base; set cntrparam levels 3; set palette rgbformulae 33,13,10; set pm3d map; set terminal wxt size 900,900; set key at 25,15 title 'Wave Function Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1); pause 0.00001; t=t+$speed;}"
        elif [ $plotOpt == 2 ]; then
            gnuplot -e "t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set zrange [0:0.13]; set xlabel 'Position x1'; set ylabel 'Position x2'; unset hidden3d; set pm3d at bs; set palette rgbformulae 33,13,10; set contour; unset surface;  set key title 'WaveFunction Probability Density'; set terminal wxt size 900,900; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1) w l; pause 0.00001; t=t+$speed;}"
        fi

        else
        
            gnuplot -e "t=0; tmax=$tIt; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position x1'; set ylabel 'Position x2'; set palette rgbformulae 33,13,10; set pm3d map; set terminal wxt size 900,900; set key at -4,10 title 'Wave Function Probability Density'; while(t<tmax){ splot 'DATA_dataToPlot_2D_CN.txt' index t title columnheader(1), 'DATA_trajectoriesToPlot_2D_CN.txt' index t title 'Trajectories' w points lt 5 pt 7 ps 1 lc rgb 'black' ; pause 0.00001; t=t+$speed;}"
        
        
        fi
    fi
    
	;;
	
	11)
	psiIni=""
    potential=""
    qDivs=""
    qmin=""
    qmax=""
    dt=1.0
    tIt=1
    plotOpt=1
    numTrajs=0
    dimNum=1
    mass=""
    speedUp=1
    potentialPlotFineness=0.02
    hbar=1.0
    substituteFor_i=""
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce n, the number of spatial dimensions of the problem"
    read dimNum
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space variables q[0], q[1],...,q[n] please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of q[0],..,q[n]"
    read potential
    echo "Introduce the mass for each dimension separated by commas m[0], m[1], ... m[n]"
    read mass
    echo "Introduce the number of space divisions for the grid in each dimension  separated by commas q[0], q[1], ... q[n]"
    read qDivs
    echo "Introduce the lower grid bounds qmin[0], qmin[1], ... qmin[n]"
    read qmin
    echo "Introduce the upper grid bounds qmax[0], qmax[1], ... qmax[n]"
    read qmax
    echo "Introduce the time iteration number "
    read tIt
    echo "Introduce the time step"
    read dt
    echo "Introduce the number of trajectories to compute"
    read numTrajs
    echo "Introduce hbar"
    read hbar
    echo "Introduce the desired substitution for the present i in the Schrodinger Equation for the dissipative time evolution (for ex. use 1 for a full imaginary time evolution, use J-(J -1)*eps with eps in (0,1) for a gradual linear interpolation between dissipative and standard time evolution...)"
    read substituteFor_i 
    echo "Introduce the fineness with which you wish to plot the potential energy field and the tensor product (0, 0.5). The lower the better, a good deal time-quality might be around 0.03"
    read potentialPlotFineness
    echo "The animation can be accelerated if you wish. Enter a speedUp factor (1/2/3/4/5/...10)"
    read speedUp
    echo ""
    gif=L
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " "$psiIni" 
    echo " V(x) = " $potential 
    echo " qDivs = " $qDivs  
    echo " qmin = " $qmin
    echo " qmax = " $qmax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo ""
    
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    if [[ $substituteFor_i != *"return"* ]]; then
        substituteFor_i="return ${substituteFor_i}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    
    ./EXE_codeFileGenerator_nD_XO_ZERO_CN_ABC_tDEP $dimNum "$psiIni" "$potential" "$mass" "$qDivs" "$qmin" "$qmax" $dt $tIt $option $numTrajs $potentialPlotFineness $hbar "$substituteFor_i"

    g++ -Wall -O CODE_simulator_nD_XO_ZERO_CN_tDEP_ImTimeEv.cpp -o EXE_simulator_nD_XO_ZERO_CN_tDEP_ImTimeEv -lm
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_nD_XO_ZERO_CN_tDEP_ImTimeEv
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."
    
    x1min=$(echo $qmin| cut -d',' -f 1)
    x2min=$(echo $qmin| cut -d',' -f 2)
    x1max=$(echo $qmax| cut -d',' -f 1)
    x2max=$(echo $qmax| cut -d',' -f 2)
    divsx1=$(echo $qDivs| cut -d',' -f 1)
    divsx2=$(echo $qDivs| cut -d',' -f 2)
    
    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  
        
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel XO_TensorProduct_$currentDateTime dim$dimNum qDivs($qDivs) qmin($qmin) qmax($qmax) numTrajs($numTrajs) dt$dt tIt$tIt mass($mass) hbar$hbar potentialPlotFineness$potentialPlotFineness speed$speedUp.gif'; t=0; countTrajs=0; tmax=2*$tIt*$numTrajs; dx1=($x1max-($x1min))/$divsx1.0; dx2=($x2max-($x2min))/$divsx2.0; array posx1[$divsx1+1]; array posx2[$divsx2+1]; do for [i=1:($divsx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($divsx2+1)] { posx2[i] = $x2min + i*dx2 }; while(t<tmax){  set multiplot; set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.025]; set xlabel 'Position q1'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index t using (posx1[\$0+1]):1 title 'q1 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 1:3 title 'q1 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$x2min:$x2max]; set yrange [0:0.025]; set xlabel 'Position q2'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t+1) using (posx2[\$0+1]):1 title 'q2 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 2:3 title 'q2 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; if ((t-(2*$tIt*countTrajs)) >= 2*$tIt){ countTrajs=countTrajs+1; }; set pm3d map; set palette rgbformulae 33,13,10; splot 'DATA_potentialToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index ((t-(2*$tIt*countTrajs))/2) title 'Potential Energy Map', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::1::(t/2+1) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';
        set origin 0.49, 0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; set palette rgbformulae 33,13,10; splot 'DATA_tensorProduct_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t/2);
        t=t+2*$speedUp; unset multiplot; }"
        echo "Done!"
    else
        gnuplot -e "t=0; countTrajs=0; tmax=2*$tIt*$numTrajs; dx1=($x1max-($x1min))/$divsx1.0; dx2=($x2max-($x2min))/$divsx2.0;  set terminal wxt size 1850,970; set multiplot; set pm3d map; array posx1[$divsx1+1]; array posx2[$divsx2+1]; do for [i=1:($divsx1+1)] { posx1[i] = $x1min + i*dx1 }; do for [i=1:($divsx2+1)] { posx2[i] = $x2min + i*dx2 }; while(t<tmax){ 
        set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [0:0.025]; set xlabel 'Position q1'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index t using (posx1[\$0+1]):1 title 'q1 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 1:3 title 'q1 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$x2min:$x2max]; set yrange [0:0.025]; set xlabel 'Position q2'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t+1) using (posx2[\$0+1]):1 title 'q2 Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::(t/2)::(t/2+1) using 2:3 title 'q2 component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; if ((t-(2*$tIt*countTrajs)) >= 2*$tIt){ countTrajs=countTrajs+1; }; set palette rgbformulae 33,13,10; splot 'DATA_potentialToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' index ((t-(2*$tIt*countTrajs))/2) title 'Potential Energy Map', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_nD_XO_ZERO_CN_ABC_tDEP.txt' every ::1::(t/2+1) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';
        set origin 0.49, 0; set size 0.5, 0.5; clear; set xrange [$x1min:$x1max]; set yrange [$x2min:$x2max]; set xlabel 'Position q1'; set ylabel 'Position q2'; set palette rgbformulae 33,13,10; splot 'DATA_tensorProduct_nD_XO_ZERO_CN_ABC_tDEP.txt' index (t/2);
        t=t+2*$speedUp; } unset multiplot; pause 3;"
            
    	#unset view; set pm3d;
    fi
    
	;;
	
	12)
	psiIni=""
    potential=""
    xDivs=0
    ydivs=0
    xmin=0.0
    xmax=0.0
    ymin=0.0
    ymax=0.0
    dt=0.01
    tIt=1
    numTrajs=0
    dimNum=1
    massx=1.0
    massy=1.0
    speedUp=1.0
    potentialPlotFineness=0.013
    hbar=1.0
    
    echo ""
    echo " (1) DATA INPUT"
    echo "Introduce in C++ syntax the initial WaveFunction as a function of space variables x, y please"
    read psiIni
    echo "Now introduce the Potential Energy field as a function of x, y"
    read potential
    echo "Introduce the mass for each dimension mx my"
    read massx massy
    echo "Introduce the number of space divisions for the grid in each dimension  xDivs yDivs"
    read xDivs yDivs
    echo "Introduce the lower grid bounds xmin ymin"
    read xmin ymin
    echo "Introduce the upper grid bounds xmax ymax"
    read xmax ymax
    echo "Introduce the time iteration number "
    read tIt
    echo "Introduce the time step"
    read dt
    echo "Introduce the number of trajectories to compute"
    read numTrajs
    echo "Introduce hbar"
    read hbar
    echo "Write in C++ syntax the shape of the eigenstates of sections in x as a function of x, y and j (the energy level)"
    read trash
    read eigenstatesForSectionsInx
    echo "Write the derivative in y of those eigenstates"
    read trash
    read diffyEigenstatesForSectionsInx
    echo "Write an expression for the second derivative in y of those eigenstates"
    read trash
    read diffyyEigenstatesForSectionsInx
    echo "Write in C++ syntax the shape of the eigenstates of sections in x as a function of x, y and j (the energy level)"
    read trash
    read eigenstatesForSectionsIny
    echo "Write the derivative in x of those eigenstates"
    read trash
    read diffxEigenstatesForSectionsIny
    echo "Write an expression for the second derivative in x of those eigenstates"
    read trash
    read diffxxEigenstatesForSectionsIny
    echo "Type the maximum j energy level to use in teh calculations for x and for y: xjmax yjmax"
    read xjmax yjmax
    echo "Give the coefficients a, b you want to use in each of the dimensions as: a_x b_x a_y b_y"
    read trash
    read b_y
    echo "Introduce the ratio of output - if 5 then only one time iteration per 5 will be outputed to the data file -> calculations and plot are faster if less data is asked to plot"
    read outputEvery
    echo "Introduce the fineness with which you wish to plot the potential energy field  (0, 0.5). The lower the better, a good deal time-quality might be around 0.03 (faster if 0.013)"
    read potentialPlotFineness
    echo "The animation can be accelerated if you wish. Enter a speedUp factor (1/2/3/4/5/...10)"
    read speedUp
    echo ""
    gif=L
    echo "If you wish to output the animation as a GIF write G, if you want it to be animated live write L"
    read gif
    
    if [[ $gif == *"G"* ]]; then
        expLabel=""
        echo "Introduce a name for the experiment - the gif file will include this label in its name besides the simulation parameters"
        read expLabel
    fi
    
    echo "-----------------"
    echo " Data received:"
    echo ""
    echo " Initial WF0(x) = " "$psiIni"
    echo " V(x) = " "$potential" 
    echo " xDivs = " $xDivs  
    echo " xmin = " $xmin
    echo " xmax = " $xmax 
    echo " yDivs = " $yDivs  
    echo " ymin = " $ymin
    echo " ymax = " $ymax 
    echo " dt = " $dt "; Number of time iterations = " $tIt
    echo ""
    
    if [[ $psiIni != *"return"* ]]; then
        psiIni="return ${psiIni}"
    fi
    if [[ $potential != *"return"* ]]; then
        potential="return ${potential}"
    fi
    
    echo " (2) CODE COMPILATION and EXECUTION"
    echo " Generating code and compiling..."
    
    ./EXE_codeFileGenerator_2D_XO_KINADV_BornHeun "$psiIni" "$potential" $massx $massy $xDivs $yDivs $xmin $xmax $ymin $ymax $dt $tIt $numTrajs $potentialPlotFineness $hbar $outputEvery "$eigenstatesForSectionsInx" "$diffyEigenstatesForSectionsInx" "$diffyyEigenstatesForSectionsInx" "$eigenstatesForSectionsIny" "$diffxEigenstatesForSectionsIny" "$diffxxEigenstatesForSectionsIny" $xjmax $yjmax $b_y
    
    g++ -Wall -O CODE_simulator_2D_XO_KINADV_BornHeun_tINDEP.cpp -o EXE_simulator_2D_XO_KINADV_BornHeun_tINDEP
    echo " Done!"
    echo ""
    echo " Executing calculations..."
    START_TIME=$SECONDS
    ./EXE_simulator_2D_XO_KINADV_BornHeun_tINDEP
    CALCULATION_TIME=$SECONDS
    echo " Done! " $(($CALCULATION_TIME - $START_TIME)) " seconds required!"
    echo ""
    
    echo " (3) PLOT ANIMATION "
    echo " Plotting..."
    
    
    if [[ $gif == *"G"* ]]; then
        currentDateTime=`date +"%m-%d-%Y_%H:%M:%S"`
        echo " Generating Animation GIF..."  
                
        gnuplot -e "set terminal gif size 1800,900 animate delay 1; set output './ANIMATION_GIFS/$expLabel XO_2D_AdvKin_CorrelPot_$currentDateTime  a,b_x($a_x,$b_x) a,b_y($a_y,$b_y) x,yjmax($xjmax,$yjmax) x,yDivs($xDivs,$yDivs) x,ymin($xmin,$ymin) x,ymax($xmax,$ymax) nTrjs($numTrajs) dt($dt) tIt($tIt) mx,y($massx,$massy) outEvry($outputEvery).gif'; t=0; countTrajs=0; tmax=2*$tIt*$numTrajs/$outputEvery; dx1=($xmax-($xmin))/$xDivs; dx2=($ymax-($ymin))/$yDivs; array posx1[$xDivs+1]; array posx2[$yDivs+1]; do for [i=1:($xDivs+1)] { posx1[i] = $xmin + i*dx1 }; do for [i=1:($yDivs+1)] { posx2[i] = $ymin + i*dx2 }; while(t<tmax){  set multiplot; set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$xmin:$xmax]; set yrange [0:0.12]; set xlabel 'Position x'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' index t using (posx1[\$0+1]):1 title 'x Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' every ::(t/2)::(t/2+1) using 1:3 title 'x component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$ymin:$ymax]; set yrange [0:0.12]; set xlabel 'Position y'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' index (t+1) using (posx2[\$0+1]):1 title 'y Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' every ::(t/2)::(t/2+1) using 2:3 title 'y component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$xmin:$xmax]; set yrange [$ymin:$ymax]; set xlabel 'Position x'; set ylabel 'Position y'; if ((t-(2*$tIt*countTrajs)) >= 2*$tIt){ countTrajs=countTrajs+1; }; set pm3d map; set palette rgbformulae 33,13,10; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' index ((t-(2*$tIt*countTrajs))/2) title 'Potential Energy Map', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' every ::1::(t/2+1) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';
        t=t+2*$speedUp; unset multiplot; }"
        
    else
        
        
        gnuplot -e "set terminal WTX size 1800,900; t=0; countTrajs=0; tmax=2*$tIt*$numTrajs; dx1=($xmax-($xmin))/$xDivs; dx2=($ymax-($ymin))/$yDivs; array posx1[$xDivs+1]; array posx2[$yDivs+1]; do for [i=1:($xDivs+1)] { posx1[i] = $xmin + i*dx1 }; do for [i=1:($yDivs+1)] { posx2[i] = $ymin + i*dx2 }; while(t<tmax){  set multiplot; set origin 0.0, 0.0; set size 0.5, 0.5; clear; set xrange [$xmin:$xmax]; set yrange [0:0.1]; set xlabel 'Position x'; set ylabel 'Probab Density';plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' index t using (posx1[\$0+1]):1 title 'x Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' every ::(t/2)::(t/2+1) using 1:3 title 'x component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.5, 0.5; set size 0.5, 0.5; clear; set xrange [$ymin:$ymax]; set yrange [0:0.1]; set xlabel 'Position y'; set ylabel 'Probab Density'; plot 'DATA_probabilityToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' index (t+1) using (posx2[\$0+1]):1 title 'y Conditional Wave Function Probability Density' w l, 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' every ::(t/2)::(t/2+1) using 2:3 title 'y component of the Trajectory' w points lt 5 pt 6 ps 2 lc rgb 'red';
        set origin 0.0, 0.5; set size 0.5, 0.5; clear; set xrange [$xmin:$xmax]; set yrange [$ymin:$ymax]; set xlabel 'Position x'; set ylabel 'Position y'; if ((t-(2*$tIt*countTrajs)) >= 2*$tIt){ countTrajs=countTrajs+1; }; set pm3d map; set palette rgbformulae 33,13,10; splot 'DATA_potentialToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' index ((t-(2*$tIt*countTrajs))/2) title 'Potential Energy Map', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' title 'Resultant Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'black', 'DATA_trajectoriesToPlot_2D_XO_CN_KinAdv_BornHeun_tINDEP.txt' every ::1::(t/2+1) title 'Trajectories' w points lt 5 pt 4 ps 0.2 lc rgb 'red';
        t=t+2*$speedUp; unset multiplot; }"

	fi
	;;
	
	*) 
    echo ""
    echo "SIMULTATION END"
    echo "Thanks for trusting XOA engine"
    echo ""
    
    break
    ;;
    
    esac

done


