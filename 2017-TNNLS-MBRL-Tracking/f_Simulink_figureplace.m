% This function places figures in a grid on the screen.

function f_Simulink_figureplace
FigHandles=find_system(gcs,'RegExp','on','Name','Scope.*');
NumFigures = size(FigHandles,1);

temp1=1:ceil(sqrt(NumFigures));
[temp2,temp3]=ndgrid(temp1,temp1);
temp1=temp2(:);
temp2=temp3(:);
temp3=temp1.*temp2-NumFigures*ones(size(temp1));
temp3(temp3<0)=inf;
[~,I]=min(temp3);
NumRows=min(temp1(I),temp2(I));
NumColumns=max(temp1(I),temp2(I));

ScnSize = get(0,'MonitorPositions');
% If there is only one display, the figures are arranged on one half
% of the display
if size(ScnSize,1)==1
    rh=ScnSize(3)/2;
    tem=NumRows; NumRows=NumColumns; NumColumns=tem;
else
    rh = ScnSize(1,3);
end
rv = ScnSize(1,4);
for i = 1:NumRows
    for j = 1:NumColumns
        if j+NumColumns*(i-1)<=NumFigures
            set_param(FigHandles{j+NumColumns*(i-1)},'Open','on');
            shh = get(0,'ShowHiddenHandles');
            set(0,'ShowHiddenHandles','On');
            set(gcf,'OuterPosition',...
                [round((j-1)*rh/NumColumns) round((NumRows-i)*rv/NumRows)...
                round(rh/NumColumns) round(rv/NumRows)]);
            set(0,'ShowHiddenHandles',shh);
        end
    end
end