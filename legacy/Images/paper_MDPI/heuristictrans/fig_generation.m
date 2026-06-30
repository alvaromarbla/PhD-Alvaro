FIG = hgload("RevHeur_Xvst.fig");
chFIG = get(FIG,"Children");
datalog = get(chFIG,"Children");
X1 = datalog{2}(1).XData;
X2 = datalog{2}(2).XData;

%YMatrix1 =  [datalog{2}(4).YData;datalog{2}(3).YData];
%YMatrix2 =  [datalog{2}(1).YData;datalog{2}(2).YData];

Y1 = datalog{2}(1).YData;
Y2 = datalog{2}(2).YData;
