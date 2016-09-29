clc
clear all;
range1 = 'A2:AH308';
range2 = 'A2:AH308';

file1 = xlsread('C:\Users\Mohanty\Documents\MATLAB\Anova\camxx_anova\monotracking_cam1.xls.ods',1, range1);
file2 = xlsread('C:\Users\Mohanty\Documents\MATLAB\Anova\camxx_anova\monotracking_cam2.xls.ods',1, range2);
filename = 'C:\Users\Mohanty\Documents\MATLAB\Anova\testdata.xlsx';

std_dev = 0.0073;           %Standard deviation in height direction (measured from calibration)
rng(0,'twister');           %Seed for normal random distribution
n = 100;  %No. of pts for distribution
K = 10;  %Max no. of groups (datasets)
X = 20;  %Factor of std Dev for defining the span   X*std_dev  defining total deviation expected
% Headers = {'Index', 'Camera1', 'Camera2', 'Overlap'};
% ToWrite = [];
% xlswrite(filename, Headers);

%CAMERA TRANFORMATION MATRIX FOR TRIANGULATION
S_0 = [208.3000     ,    0       ,  0;
         0,  208.3000   ,      0];
S_1 = [181.1000     ,    0       ,  0;
         0,  181.8000   ,      0];
W_0 = [-1.1991;-1.2237;1.2150];
W_1 = [0.0127;0.0097;1.5813];
T_0 = [4.6392;2.9923];
T_1 = [5.6624;5.5949];

%Buffers for storing the triangulated data points
 X_first =[];
 X_second =[];
 Prob = [];

       

%Choosing 5 detections from camera 1
for i = 1:307
    
idx =4;
C1 = [];
count1 = 0;
    while ((idx<34)&& ~isnan(file1(i,idx))) 
        C1 = [C1, X*std_dev.*(randn(n,1))+file1(i,idx)];
        idx = idx + 7;
        count1 = count1+1;
    end
    
idx =4;
C2 = [];
count2 = 0;

    %Chossing 5 detections from camera 2
    while ((idx<34)&& ~isnan(file2(i,idx))) 
        C2 = [C2, X*std_dev.*(randn(n,1))+file2(i,idx)];
        idx = idx + 7;
        count2 = count2+1;
    end
    
K = count1+count2;      %Total number of points (camera 1+2) for comparison

C = [C1,C2];            %Augmented dataset 
%figure(1)
[p, tbl, stats] = anova1(C,[], 'off');          %ANOVA
% figure(2)
% title('ANOVA test for triangulation');
% set(gca, 'FontSize', 12, 'FontName', 'TimeNewRoman', 'FontWeight', 'bold')
% ylabel('Height (mm)');
% xlabel('Point number');
% %%
[c,~,~,gnames] = multcompare(stats);            %Tukey-HSD

% figure(1)
% title('Tukey-Kramer HSD test for triangulation');
% set(gca, 'FontSize', 12, 'FontName', 'TimeNewRoman', 'FontWeight', 'bold')
% xlabel('Height (mm)');
% ylabel('Point number');

% WritePoint = [];

for j = 1:(K*(K-1)/2)       %Performing contrast of means for all points in tracking of trajectory
if (c(j,6)>0.1)             %Probability of overlap in distribution
    if (c(j,1)>0 && c(j,1)<count1+1) && (c(j,2)>K-count2 && c(j,2)<K+1)         %Point1 from Camera1 ; Point2 from Camera2
       
     k1 = find(abs(file1(i,:)-stats.means(c(j,1)))<0.1);   %Idx of pixels in cam1
     k2 = find(abs(file2(i,:)-stats.means(c(j,1)))<0.1);   %Idx of pixels in cam1
    disp(['Point ', num2str(i) ,  ' : Camera 1 : ', num2str(c(j,1)), ', Camera 2 : ', num2str(c(j,2)), ', Significance : ', num2str(c(j,6)) ,' Height : ' ,num2str(stats.means(c(j,1))),' , ', num2str(stats.means(c(j,2)))]);
    ToWrite = [ ToWrite; i , c(j,1), c(j,2), c(j,6), stats.means(c(j,1)), stats.means(c(j,2)), file1(i,k1(1)-2), file1(i,k1(1)-1), file2(i,k2(1)-2), file2(i,k2(1)-1)];
    %ToWrite = [ToWrite;WritePoint];
    xlswrite(filename, ToWrite);
    
    %Triangulating the points obtained with given probability
    UV_first = [file1(i,k1(1)-2); file1(i,k1(1)-1); file2(i,k2(1)-2); file2(i,k2(1)-1)];
    Prob_next = (c(j,6));
    %UV_second = [file2(i,k2(1)-2) file2(i,k2(1)-1)]
    X_first_next = ([S_0*Rodrigues(W_0);S_1*Rodrigues(W_1)])\(UV_first - [S_0,zeros(2,3);zeros(2,3),S_1]*[T_0;0;T_1;0]);
    %X_second_next = ([S_0*Rodrigues(W_0);S_1*Rodrigues(W_1)])\(UV_second - [S_0,zeros(2,3);zeros(2,3),S_1]*[T_0;0;T_1;0]);
              
    X_first = [X_first, X_first_next];
    %X_second = [X_second, X_second_next];
    Prob = [Prob, Prob_next];
             
    break
    end
end
end
end

%Transforming coordinates from smarAct to assumed world frame
x_0 = vrrotmat2vec([0,-1,0;0,0,1;-1,0,0]);
angleAxis = [x_0(1:3)*x_0(4)]';
X_first = Rodrigues(angleAxis)*X_first;
%X_second = Rodrigues(angleAxis)*X_second;
            
%%

%Plotting the triangulated points based on probability of overlap
figure(1)
 grid on
 for i = 1:size(X_first, 2)
     xlabel('X (mm)');
     ylabel('Y (mm)');
     zlabel('Z (mm)');
     ax = gca;
    ax.XLim = [-4 4];
    ax.YLim = [-4 4];
    ax.ZLim = [-4 4];
  
    if(Prob(i)<0.5)
   
    h = scatter3(ax,  X_first(1,i), X_first(2,i), X_first(3,i), [], [1,0,0]);
    end
    
    if(Prob(i)>0.5 && Prob(i)<0.9)
    h = scatter3(ax,  X_first(1,i), X_first(2,i), X_first(3,i), [], [0,1,0]);
    end
    
    if(Prob(i)>0.9)
    h = scatter3(ax,  X_first(1,i), X_first(2,i), X_first(3,i), [], [0,0,1]);
    end
    hold on
    %h = scatter3(ax,  X_second(1,i), X_second(2,i), X_second(3,i), [], [0,0,1]);

  
    title('Trajectory of ABFs in workspace');
    set(gca, 'FontSize', 12, 'FontName', 'TimeNewRoman', 'FontWeight', 'bold')
     box on
     drawnow
     hold on
     

 end
