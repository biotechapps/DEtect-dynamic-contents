clc;
clear all;
prompt = 'Please enter the Webpage start time(sec) ';
count = input(prompt);

prompt = 'Enter the Webpage end-session (sec) ';
Time = input(prompt);
 d = count:Time;

prompt = 'Enter the Webpage Name or URL: ';
I = input(prompt,'s');

format compact 
 figure('Menubar','none','Name','Eye Movement prediction', ... 
         'NumberTitle','off'); 
I = imread(I);
  

PP12 = dataset('xlsfile','PP12.xlsx');
x = PP12.MappedFixationPointX(d);
y = PP12.MappedFixationPointY(d);
z = PP12.FixationIndex(d);

SCR = PP12.eda(d);
ST = PP12.ST(d);
MappedFX = PP12.MappedFixationPointX(d);
MappedFY = PP12.MappedFixationPointY(d);
Pupil_changes  = PP12.Pupil_changes(d);

data3 = [SCR ST MappedFY Pupil_changes];
data4 = [SCR ST MappedFX Pupil_changes];

data3 = iddata(x,data3,0.08);
data3 = misdata(data3);
MDL3 = PHYCOB(data3, 3,'StateName',{'order1','order2','order3'},'OutputName',{'MappedFX'},'InputName',{'SCR','ST','MappedFY','Pupil_changes'});

data4 = iddata(y,data4,0.08);
data4 = misdata(data4);
MDL4 = PHYCOB(data4, 3,'StateName',{'order1','order2','order3'},'OutputName',{'MappedFY'},'InputName',{'SCR','ST','MappedFX','Pupil_changes'});
size(MDL3);

K = 4;
[ypp,x0,syspp] = predict(MDL3,data3,K);
[yppp,x0,sysppp] = predict(MDL4,data4,K);

opt = simOptions;
opt.InitialCondition = x0;
ys = sim(syspp,[data3.OutputData,data3.InputData],opt);

pred2 = round(ypp.OutputData);
pred22 = round(yppp.OutputData);

 [peak, locs] = findpeaks(SCR,'minpeakdistance',15);
            m = length(locs);
               emotn1 = strmatch('Stressed',PP12.Affectstate(locs(1:m)));
             emotn2 = strmatch('Relaxed',PP12.Affectstate(locs(1:m)));    
            if size(emotn2) == [0,1];
    emotn2 = 1
            else
    emotn2 = emotn2
            end
         
            X = PP12.MappedFixationPointX(locs(emotn1));
            Y = PP12.MappedFixationPointY(locs(emotn1));
            
           XX =PP12.MappedFixationPointX(locs(emotn2));
          YY = PP12.MappedFixationPointY(locs(emotn2));                          
               


xi = x(1);
yi = y(1);
zi = z(1);
xImg = linspace(min(x), max(x), size(I, 2));
yImg = linspace(min(y), max(y), size(I, 1));
image(xImg, yImg, I, 'CDataMapping', 'scaled');
hold on;

h = plot(xi,yi,'mo', ...
         xi, yi, 'r-','MarkerEdgeColor','w',...
	 'MarkerFaceColor','m',...
         'MarkerSize',23,'color','m');

ylabel('MappedFixationPointY');
xlabel('MappedFixationPointX');
axis off

for i=1:numel(x);
htext = text(x(i),y(i),[' ', num2str(z(i)), ''],...
			 'Color',[0 0 0], 'FontSize',8, ...
                     'HorizontalAlignment','center', 'color','w','VerticalAlignment','middle');
       title(['elapse time ' num2str(i)])

set(htext, 'Position',[x(i) y(i)])
set(h,'XData',x(1:i),'YData',y(1:i))              
    drawnow; pause(0.5);
end   

hold on
D=rand(length(pred2),length(pred2));
for i=1:length(D(:,:));
   scatter(pred2,pred22,(700*D(:,i)),'filled','MarkerFaceAlpha',.3);
        drawnow
    pause(0.01)
end

hold on;
             for i=1:length(D(:,1));
    s = scatter(X, Y,'filled','MarkerFaceColor','r','SizeData',30000);
    alpha(s,0.01)
    
    scatter(XX,YY,'filled','SizeData',30000,'MarkerFaceAlpha',0.01);
    drawnow
    pause(0.01)
             end
             
             % MDL_2nd_order = idpoly(MDL4);
      MDL4 = frd(MDL4,SCR);
%bode(mdl,'r-');
MDL_2nd_order = fitfrd(MDL4,2);
MDL_1st_order = fitfrd(MDL4,1);
MDL_3rd_order = fitfrd(MDL4,3);
bode(MDL_2nd_order,'r',MDL_1st_order,'k:',MDL_3rd_order,'b-.');
legend('Show');
%print('yahoo','-dpdf','-fillpage')
%Gmot = zpk([],[-1,-1],1);
%Cmot = tunablePID('Cmot','PID');
%controlSystemTuner