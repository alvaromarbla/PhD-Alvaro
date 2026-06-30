FIG = hgload("EpervstOC50.fig");
chFIG = get(FIG,"Children");
datalog = get(chFIG,"Children");
X1 = datalog{2}(4).XData;
X2 = [datalog{2}(1).XData,datalog{2}(2).XData,datalog{2}(3).XData];

YMatrix1 = [datalog{2}(6).YData;datalog{2}(4).YData;datalog{2}(5).YData];
YMatrix2 = [datalog{2}(3).YData;datalog{2}(1).YData;datalog{2}(2).YData];

