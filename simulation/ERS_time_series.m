function acqdates = ERS_time_series(infolevel,N)
%
% acqdates = ERS_time_series(infolevel,[N])
%
% simulate a distribution of data using ERS characteristics 
% INPUT:
% N (optional): number of simulated ERS acquisition dates
% infolevel=0 : no extra output, 
% infolevel=1 : no graphical output, only table 
% infolevel=2 : no graphical output, more information on screen, and table
% infolevel=3 : also graphical information output
% infolevel=4 : no output at all, show only ERS data acquisition phases
%
% OUTPUT:
%   acqdates = [acquisition dates , sat_ID]
%                in nr of days since  1-Jan-0000
%
% Example  acqdates = ERS_time_series(3,20)
%  

% Ramon Hanssen nov 2002

 
if nargin == 0, help ERS_time_series;return;end
if nargin == 1, 
  if infolevel == 4, 
     simulation = 0; 
     fprintf(1,'Showing only ERS acquisition phases\n')
  else  
      fprintf(1,'\n>>>Give number of image acquisitions N to simulate<<<\n\n');
      help ERS_time_series;
      return
  end 
  acqdates=[];
end
if nargin==2,
  simulation = 1; 
end

%infolevel=3,N=80,simulation=1

%% from http://www.deos.tudelft.nl/ers/phases/
%rep.freq.
%E1    3   Phase A: Commissioning Phase: 25-Jul-91/10-Dec-91 
%E1    3   Phase B: Ice Phase: 28-Dec-91/30-Mar-92 
%E1   35   Phase R: Roll Tilt Mode Campaign: 4-Apr-92/13-Apr-92 
%E1   35   Phase C: Multidisciplinary Phase: 14-Apr-92/21-Dec-93 
%E1    3   Phase D: Second Ice Phase: 23-Dec-93/10-Apr-94 
%E1  168   Phase E: First Geodetic Phase:10-Apr-94/28-Sep-94 
%E1  168   Phase F: Second Geodetic Phase: 28-Sep-94/21-Mar-95 
%E1   35   Phase G: Tandem Phase: 21-Mar-95/5-Jun-96 
%E2   35   Phase A: Multi-disciplinary Phase: 29-Apr-96/now 
%
%E1   35   Phase R: datenum('4-Apr-92');datenum('13-Apr-92')
%E1   35   Phase C: datenum('14-Apr-92');datenum('21-Dec-93') 
%E1   35   Phase G: datenum('21-Mar-95');datenum('5-Jun-96') 
%E2   35   Phase A: datenum('29-Apr-96');datenum(date)


ERS1_start = datenum('25-Jul-91');
ERS2_start = datenum('21-Mar-95');
%ERS2_start = datenum('29-Apr-96');
TODAY      = datenum(date);
ERS1_R     = datenum('4-Apr-92'):datenum('13-Apr-92');
ERS1_C     = datenum('14-Apr-92'):datenum('21-Dec-93') ;
ERS1_G     = datenum('21-Mar-95'):datenum('5-Jun-96');
ERS2_tandem= datenum('21-Mar-95'):datenum('5-Jun-96');
ERS2_A     = ERS2_start:TODAY;


Lifetime   = [ERS1_start:TODAY];
ERS1_on     = [ERS1_R,ERS1_C,ERS1_G];
ERS2_on     = [ERS2_tandem,ERS2_A];
[c,e1a,e1b]  = intersect(Lifetime,ERS1_on);
[c,e2a,e2b]  = intersect(Lifetime,ERS2_on);
Life         = Lifetime-min(Lifetime);  % days from start ERS-1
operation1  = zeros(size(Lifetime));
operation2  = zeros(size(Lifetime));
operation1(e1a) = 1;
operation2(e2a) = 0.95;

if ( infolevel== 3 | infolevel== 4)
   figure(1);plot(Life,operation1,'b','linewidth',[4]);
   set(gca,'ylim',[0 2],'ytick',[0 1])
   xlabel('Days from start ERS-1');
   ylabel('Operation yes (1) no (0)');
   hold on;
   plot(Life,operation2,'m','linewidth',[4]);
   for i=1991:2004,
     y=datenum([i,01,01,0,0,0])-min(Lifetime);
     h=line([y y],[0 2],'linestyle',':');
     h=text(y,1.5,num2str(i),'rotation',[90]);
   end
   hold off
   legend('ERS-1','ERS-2')
end

%%% Simulation part
if simulation == 1
   % N is the number of simulations
   RI = 35;  % repeat interval
   % all on days
   ERSboth_on = unique([ERS1_on,ERS2_on]);

   % create a random reference date
    qqq=ceil(rand*(length(ERSboth_on)));
   refdate  = ERSboth_on(qqq);

   % select all possible acquisition days wrt the refdate
   E1_tot      = [refdate:RI:max(ERS1_on),refdate:-RI:ERS1_start];
   E2_tot      = [refdate+1:RI:max(ERS2_on),refdate+1:-RI:ERS2_start];
   E1E2_tot    = sort ([E1_tot,E2_tot]);
   E1E2_tot_random_order  = E1E2_tot(randperm(length(E1E2_tot)));
   E1E2_select = sort(E1E2_tot_random_order(1:N));
   satellite_ID= zeros(size( E1E2_select));
   satellite_ID(ismember(E1E2_select,E1_tot))=1;
   satellite_ID(find(satellite_ID==0))=2;
   acqdates=([E1E2_select',satellite_ID']);
   if ( infolevel == 2 | infolevel == 3 )
     fprintf(1,'SAT\tDATE\tDAYNR\n')
     for jj=1:length(acqdates(:,1)),
       fprintf(1,'ERS-%i \t %s \t%i\n',acqdates(jj,2),datestr(datevec(acqdates(jj,1))),acqdates(jj,1))
     end
   end
   if ( infolevel== 3)
      figure(1);hold on
      %plot(E1_tot-min(Lifetime),ones(size(E1_tot)),'bo')
      %plot(E2_tot-min(Lifetime),ones(size(E2_tot)),'mo')
      plot(E1E2_select-min(Lifetime),ones(size(E1E2_select)),...
          'ro','markersize',10)
      %plot(E1E2_tot-min(Lifetime),ones(size(E1E2_tot)),'o')
      hold  off
   end
end
