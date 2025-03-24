unit FormMain;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, TAGraph, TASeries, Forms, Controls, Graphics,
  Dialogs, StdCtrls, ComCtrls, math, types, TACustomSeries, TADrawUtils,
  TAMultiSeries, TATransformations;

type

  { TForm1 }

  TForm1 = class(TForm)
    ButtonStop: TButton;
    ButtonRun: TButton;
    BrightnessAreaSeries: TAreaSeries;
    ChartHistorgram: TChart;
    TrapDepthLineSeries: TLineSeries;
    TrapsAreaSeries: TAreaSeries;
    OFFBarSeries: TBarSeries;
    ONBarSeries: TBarSeries;
    HistogramBarSeries: TBarSeries;
    ChartLOGLOG: TChart;
    TrapsLineSeries: TLineSeries;
    LabelKT: TLabel;
    kTTrackBar: TTrackBar;
    TrionsCountLineSeries: TLineSeries;
    ECountLineSeries: TLineSeries;
    ChartTrace: TChart;
    procedure ButtonRunClick(Sender: TObject);
    procedure ButtonStopClick(Sender: TObject);
    procedure kTTrackBarChange(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end;

  TTrap=Record
    x,y:integer;
    isActive:boolean;
    Energy:real;
  end;

  TSite=Record
    Energy:real;
    isTrap:boolean;
  end;





const
  NTraps=  20;
  Size =   100;           // size of the grid on which traps/adatoms travel
  TrapFraction = 1/Ntraps;   // if trap probability is 50%, then out of 10 traps, 5 would always be active, making it permanently dark
  //TrapDepth = 0.5;     // eV. Smaller values makes it more flickery (intermediate intensities in histogram)
  //barrierFraction = 0.35;// Marcus, Scholes PNAS. deltaG = 130 meV, reorg energy = 200 meV (i.e. jump forward = 200 meV, jump backward = 70 meV = 0.35 * 200)
  //deltaG = 1;          // eV. Smaller values makes it more flickery (intermediate intensities in histogram)
  //lambda_forward  = 1; //200/130;         // reorganization energy, Marcus, Scholes PNAS. deltaG = 130 meV, reorg energy = 200 meV (i.e. jump forward = 200 meV, jump backward = 70 meV = 0.35 * 200)
  //lambda_backward = (200 - 130)/130;      // reorganization energy, Marcus, Scholes PNAS. deltaG = 130 meV, reorg energy = 200 meV (i.e. jump forward = 200 meV, jump backward = 70 meV = 0.35 * 200)
  FWHM = 0.5;              // width of the Gaussian for trap energy distribution. Scholes PNAS, sigma 50 meV i.e. width 100 meV, vs 130 meV depth
  DEPTH_SCALING = 0.05;   // 2x higher barrier decreases the jumping chance 7.4x, before was 1.5 let me make it 3

  //kT = 300 * 0.027/300;  // eV

  bin =    10e-3;        // 10 ms
  subbins = 350;         // resolution up to 300 photons per bin
  timestep=bin/subbins;

  k_rad =      1/20e-9;  // tau = 15-30 ns
  k_trapping = 1/0.1e-9; // trapping competes with radiative decay. Shell effect is very dramatic

  k_abs =   250/bin;     // observed brightness is usually ~200 photons/bin, in reality, due to flickering, number of absorbed photons are likely much higher
  k_abs_bg = 3/bin;
  threshold = 60;

  temp = 300;

  version = 1;





var
  Form1: TForm1;

  i,j,k,trap_i, Ntrapped: integer;
  canStop: boolean;
  time: integer;
  traps: array[1..NTraps]of TTrap;
  sites: array[1..NTraps, 1..Size, 1..Size]of TSite;
  rate: array[1..10]of real;
  histogram: array[1..350] of integer;
  ONevents, OFFevents: array[1..50000] of integer;   // up to 100 seconds duration
  kT: real;
  SelectedTraps: double;
  MinTrap: integer;
  MaxTrap: integer;
  //TrapFraction: double;





implementation

{$R *.lfm}

{ TForm1 }




procedure TForm1.ButtonRunClick(Sender: TObject);
var
  photonsAbsorbed, photonsEmitted: integer;
  TrapsActive: integer;
  ONduration, OFFduration: integer;
  isON, wasON: boolean;
  sum, gamma_plus, gamma_minus: real;
  prev_i, prev_j, next_i, next_j: integer;
  event, prob_position: real;
  event_index: integer;
  k_nonrad: real;
  logFile_Photons: TextFile;
  logFile_Histogram: TextFile;
  logFile_Sites: TextFile;
  logFileName_Photons: string;
  logFileName_Histogram: string;
  logFileName_Sites: string;
  FWHMnew: Double;

begin
  Randomize;
  canStop:=false;
  time:=0;
  photonsAbsorbed:=0;
  photonsEmitted:=0;

  kT:=kTTrackBar.Position/temp*0.027;
  FWHMnew:=kT*FWHM;


  ONduration:=0;
  OFFduration:=0;
  wasON:=True;
  isON:=True;

  BrightnessAreaSeries.Clear;
  HistogramBarSeries.Clear;
  for i:=0 to high(histogram) do
    histogram[i]:=0;
  ONBarSeries.Clear;
  for i:=0 to high(ONevents) do
    ONevents[i]:=0;
  OFFBarSeries.Clear;
  for i:=0 to high(OFFevents) do
    OFFevents[i]:=0;

  //logFileName_Photons := 'Version_blinking_2025_20_300k.txt';
  logFileName_Photons := Format('Version_blinking_2025_ogFWHM_%d_%dk_depths_%.0f_fwhm_%.0f_version_%d_NTraps_%d.txt',
                                [temp div 10, temp, DEPTH_SCALING * 1000, FWHM * 1000, version, NTraps]);


  //initialize the grid of sites ===========================================================================================================
  for trap_i:=1 to NTraps do   // separate landscape of different depth for each trap!
  for i:=1 to Size do
  for j:=1 to Size do
  begin
    sites[trap_i, i, j].Energy := (FWHM * DEPTH_SCALING * trap_i / (2 * Sqrt(2 * ln(2)))) * Sqrt(-2 * ln(Random)) * cos(2 * Pi * Random); // Box-Meuler method to recreate Gaussian

    //MinTrap := NTraps; // Minimum trap
    //MaxTrap := NTraps; // Maximum trap

    SelectedTraps := Random(MaxTrap - MinTrap + 1) + MinTrap; // Random value in range [MinTrap, MaxTrap]

    //SelectedTraps := Random(Ntraps) + 1; //randomize for each site
    //TrapFraction := 1 / SelectedTraps;
    if random <= TrapFraction then
      sites[trap_i, i,j].isTrap:=True
    else
      sites[trap_i, i,j].isTrap:=False;
  end;



  //initialize trap locations on the sites grid
  for trap_i:=1 to NTraps do
  begin
    traps[trap_i].x:=random(Size)+1;  // random is from 0 to strictly <Size
    traps[trap_i].y:=random(Size)+1;

    traps[trap_i].isActive:=sites[trap_i, traps[trap_i].x, traps[trap_i].y].isTrap;
  end;









  repeat   //==============================================================================================================================
    inc(time);


    //move the adatom traps over the surface sites grid   ---------------------------------------------------------------------------------
    for trap_i:=1 to NTraps do
    begin
        //current trap location
        i:=traps[trap_i].x;
        j:=traps[trap_i].y;

        //periodic boundaries (for the upcoming step)
        prev_i:=i-1; if prev_i<1    then prev_i:=Size;
        next_i:=i+1; if next_i>Size then next_i:=1;

        prev_j:=j-1; if prev_j<1    then prev_j:=Size;
        next_j:=j+1; if next_j>Size then next_j:=1;




        // Probabilities of jumping in each direction:
        // (we don't restrict jumping into populated site, i.e. traps are independent or not floating o nthe same surface)

        // barrier to jump left                                                                                      //x-1
        if((sites[trap_i, prev_i,j].Energy - sites[trap_i, i,j].Energy)>0)then
          rate[1]:=1/timestep/4*exp(-(sites[trap_i, prev_i,j].Energy -sites[trap_i, i,j].Energy)/kT) //jump uphill, damp it exp
        else
          rate[1]:=1/timestep/4; // jump downhill, don't pre-set any barrier, they will appear just due to grid
          // divided by 4 is because we have 4 directions, and when we'll sum them up, we shouldn't overflow
          // usually Monte-Carlo is done with a constant rate downhill (no barrier), in our case it leads to a lot of flickering, which we actually do need
          // (old, WITH a barrier: //*exp(+(sites[trap_i, i,prev_j].Energy -sites[trap_i, i,j].Energy)*lambda_backward/kT);)

        //barrier to jump right                                                                                       //x+1
        if((sites[trap_i, next_i,j].Energy - sites[trap_i, i,j].Energy)>0)then
          rate[2]:=1/timestep/4*exp(-(sites[trap_i, next_i,j].Energy -sites[trap_i, i,j].Energy)/kT)
        else
          rate[2]:=1/timestep/4;

        //barrier to jump down                                                                                         //y-1
        if((sites[trap_i, i,prev_j].Energy - sites[trap_i, i,j].Energy)>0)then
          rate[3]:=1/timestep/4*exp(-(sites[trap_i, i,prev_j].Energy -sites[trap_i, i,j].Energy)/kT)
        else
          rate[3]:=1/timestep/4;

        //barrier to jump up                                                                                           //y+1
        if((sites[trap_i, i,next_j].Energy - sites[trap_i, i,j].Energy)>0)then
          rate[4]:=1/timestep/4*exp(-(sites[trap_i, i,next_j].Energy -sites[trap_i, i,j].Energy)/kT)
        else
          rate[4]:=1/timestep/4;

        //nothing happens
        rate[5]:=1e5;




        //throw the dice for the trap hop
        event:=random;
        prob_position:=0;
        for event_index:=1 to 5 do   //2x + 2y + 1nothing
        begin
          prob_position:= prob_position + timestep*rate[event_index];
          if (event<prob_position) then break;
        end;



         //what to do for the chosen event
         case event_index of // if nothing happend, last event_index would be used
         1:begin
             traps[trap_i].x:=prev_i;
             traps[trap_i].isActive:=sites[trap_i, traps[trap_i].x, traps[trap_i].y].isTrap;
             traps[trap_i].Energy:=sites[trap_i, traps[trap_i].x, traps[trap_i].y].Energy;
           end;
         2:begin
             traps[trap_i].x:=next_i;
             traps[trap_i].isActive:=sites[trap_i, traps[trap_i].x, traps[trap_i].y].isTrap;
             traps[trap_i].Energy:=sites[trap_i, traps[trap_i].x, traps[trap_i].y].Energy;
           end;

         3:begin
             traps[trap_i].y:=prev_j;
             traps[trap_i].isActive:=sites[trap_i, traps[trap_i].x, traps[trap_i].y].isTrap;
             traps[trap_i].Energy:=sites[trap_i, traps[trap_i].x, traps[trap_i].y].Energy;
           end;
         4:begin
             traps[trap_i].y:=next_j;
             traps[trap_i].isActive:=sites[trap_i, traps[trap_i].x, traps[trap_i].y].isTrap;
             traps[trap_i].Energy:=sites[trap_i, traps[trap_i].x, traps[trap_i].y].Energy;
           end;
         //5://default, do nothing
         end;



    end; //cycle for all traps















    //Sandor's model, count how many active traps right now  -------------------------------------------------------------------
    //sum:=0;
    //for i:=1 to NTraps do
    //  sum:=sum+(Ord(traps[i].isActive)-0.5); // numbers from previous step. they are not changing on the fly








    //activate and count traps  -------------------------------------------------------------------------------------------------
    TrapsActive:=0;
    k_nonrad:=0;
    for trap_i:=1 to NTraps do
    begin
       if traps[trap_i].isActive then
       begin
          inc(TrapsActive); // counter to plot in real time
          k_nonrad:=k_nonrad + k_trapping
       end;
    end;
    //TrapsLineSeries.AddXY(time/100/100, -TrapsActive*1);
    //TrapsAreaSeries.AddXY(time/100/100, TrapsActive*1); // this is a correct (but very slow) way of plotting trap flips









    //absorb light  ----------------------------------------------------------------------------------------------------------------
    if random <= k_abs * timestep then    // i.e. 250/300. Do not make it 300/300, let it have some realistic noise
       photonsAbsorbed:=1
    else
       photonsAbsorbed:=0;


    // emit light vs. get trapped (with PLQY probability)
    if random <= k_rad / (k_rad + k_nonrad) then
    begin
       photonsEmitted:=photonsEmitted + photonsAbsorbed; //not every timestep absorbs a photon!
       //photonsAbsorbed:=0;  // shoudl we do dec() instead?
    end

    else //photon is trapped. no NEW photons emitted, i.e. total photonsEmitted per bin remain unchanged
    begin
       // no photon emitted, but we don't reset it yet, because it is just 1 timestep, but the bin hasn't finished yet
       //photonsAbsorbed:=0;
    end;


    //background photons ~5/bin
    if random <= k_abs_bg * timestep then
       inc(photonsEmitted);







    AssignFile(logFile_Photons, logFileName_Photons);
    {$I-} // Turn off automatic error checking
    if FileExists(logFileName_Photons) then
        Append(logFile_Photons)
    else
        Rewrite(logFile_Photons);
    {$I+} // Turn automatic error checking back on

    if IOResult <> 0 then
        ShowMessage('Error opening or creating the file.');


    // plot results at the end of the bin (but in reality there are 100 trap events within the bin!)
    if (time mod subbins) =0 then       // i.e. every bin add a point
    begin
        BrightnessAreaSeries.AddXY(time*timestep, photonsEmitted);
        TrapsAreaSeries.AddXY(time*timestep, TrapsActive*1);  // plotting at the end of the bin makes it 100x faster but loses trap flickering within the bin
        TrapDepthLineSeries.AddXY(time*timestep, traps[1].Energy*100 + 250); // current site depth fluctuations
        Writeln(logFile_Photons, FloatToStr(time*timestep) + ' ' + FloatToStr(photonsEmitted));

        inc(histogram[photonsEmitted]);
        HistogramBarSeries.AddXY(photonsEmitted, histogram[photonsEmitted]); // ,'',time);



        if photonsEmitted > threshold then
        //if (25 < photonsEmitted) and (photonsEmitted < 65) then    // THIS DOES NOT WORK PROPERLY!!! gives long OFF states when there are no OFF states!
          isON:=True
        else
          isON:=False;


        if isON = wasON then  // event continues
          begin
            if isON then
               inc(ONduration)
            else
               inc(OFFduration);
          end

        else   // switch happened
          begin
            if wasON then
            begin
               inc(ONduration);  //if 2 switches in a row happen, we still should get duration 1
               inc(ONevents[ONduration]);
               //dec(ONevents[ONduration-1]);  don't need this because switch has actually completed
               ONBarSeries.AddXY(log10(ONduration), log10(ONevents[ONduration])+0.01);
            end

            else //wasOFF
            begin
               inc(OFFduration);
               inc(OFFevents[OFFduration]);
               //dec(OFFevents[OFFduration-1]);  don't need this because switch has actually completed
               OFFBarSeries.AddXY(log10(OFFduration)+0.001, log10(OFFevents[OFFduration])+0.01);
            end;

            // reset for the next cycle
            ONduration:= 0;
            OFFduration:=0;
          end;




        //reset for the next cyclce
        wasON:=isON;
        photonsEmitted:=0;
    end;

    CloseFile(logFile_Photons);


    if((time mod int(60*100*subbins-1))=0) then
      Application.ProcessMessages;

    if((time mod int(60*100*subbins))=0) then
    begin
      BrightnessAreaSeries.Clear;
      //TrapsLineSeries.Clear;
      TrapsAreaSeries.Clear;
      TrapDepthLineSeries.Clear;
      Form1.Caption:=IntToStr(traps[1].x)+' '+IntToStr(traps[1].y);
    end;

  until canStop or (time>3600*100*subbins);//1 hr * 100 bins * subbins photons
end;


procedure TForm1.ButtonStopClick(Sender: TObject);
begin
  canStop:=true;
end;

procedure TForm1.kTTrackBarChange(Sender: TObject);
begin
  kT:=kTTrackBar.Position/300*0.027;
end;





end.

