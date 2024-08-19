function [grid_data] = Data_reader_industrial()

filename = '..\Atlantide networks\ATL_Grid_Industrial.xls';
[~,grid_data.from] = xlsread(filename,'Linee','B2:B232');
[~,grid_data.to] = xlsread(filename,'Linee','C2:C232');
[~,grid_data.type] = xlsread(filename,'Linee','D2:D232');
grid_data.len = xlsread(filename,'Linee','F2:F232');
[~,grid_data.state_from] = xlsread(filename,'Linee','H2:H232');
[~,grid_data.state_to] = xlsread(filename,'Linee','I2:I232');



filename2 = '..\Atlantide networks\ATL_Components_Industrial_Grid.xls';
[~,grid_data.cable] = xlsread(filename2,'Linee','A2:A44');
grid_data.R = xlsread(filename2,'Linee','C2:C44');
grid_data.L = xlsread(filename2,'Linee','D2:D44');




[~,grid_data.ff] = xlsread(filename,'Nodi','F2:F101');




filename3 = '..\Atlantide networks\ATL_Profiles_Industrial_Grid.xls';
grid_data.coeff.res.day = xlsread(filename3, 'Profili Carico Giornalieri', 'C2:C97');
grid_data.coeff.ind.day = xlsread(filename3, 'Profili Carico Giornalieri', 'D2:D97');
grid_data.coeff.com.day = xlsread(filename3, 'Profili Carico Giornalieri', 'E2:E97');

grid_data.coeff.ind.week = xlsread(filename3, 'Coeff Sett Carico', 'C2:C8');
grid_data.coeff.com.week = xlsread(filename3, 'Coeff Sett Carico', 'D2:D8');
grid_data.coeff.res.week = xlsread(filename3, 'Coeff Sett Carico', 'E2:E8');

grid_data.coeff.ind.month = xlsread(filename3, 'Coeff Mens Carico', 'C2:C13');
grid_data.coeff.com.month = xlsread(filename3, 'Coeff Mens Carico', 'D2:D13');
grid_data.coeff.res.month = xlsread(filename3, 'Coeff Mens Carico', 'E2:E13');

grid_data.coeff.res.year = xlsread(filename3, 'Coeff Ann Carico', 'C2:C22');
grid_data.coeff.ind.year = xlsread(filename3, 'Coeff Ann Carico', 'D2:D22');
grid_data.coeff.com.year = xlsread(filename3, 'Coeff Ann Carico', 'E2:E22');

grid_data.coeff.pv.day = xlsread(filename3, 'Profili Gen Giornalieri', 'D2:D97');
grid_data.coeff.chp.day = xlsread(filename3, 'Profili Gen Giornalieri', 'E2:E97');
grid_data.coeff.wind.day = xlsread(filename3, 'Profili Gen Giornalieri', 'F2:F97');

grid_data.coeff.pv.week = xlsread(filename3, 'Coeff Sett Gen', 'D2:D8');
grid_data.coeff.chp.week = xlsread(filename3, 'Coeff Sett Gen', 'E2:E8');
grid_data.coeff.wind.week = xlsread(filename3, 'Coeff Sett Gen', 'F2:F8');

grid_data.coeff.pv.month = xlsread(filename3, 'Coeff Mens Gen', 'D2:D13');
grid_data.coeff.chp.month = xlsread(filename3, 'Coeff Mens Gen', 'E2:E13');
grid_data.coeff.wind.month = xlsread(filename3, 'Coeff Mens Gen', 'F2:F13');

grid_data.coeff.pv.year = xlsread(filename3, 'Coeff Ann Gen', 'D2:D22');
grid_data.coeff.chp.year = xlsread(filename3, 'Coeff Ann Gen', 'E2:E22');
grid_data.coeff.wind.year = xlsread(filename3, 'Coeff Ann Gen', 'F2:F22');



[~,grid_data.load_bus] = xlsread(filename, 'Carichi_A', 'B2:B129');
grid_data.load_P = xlsread(filename, 'Carichi_A', 'C2:C129');
grid_data.load_Q = xlsread(filename, 'Carichi_A', 'D2:D129');
[~,grid_data.load_type] = xlsread(filename, 'Carichi_A', 'F2:F129');


[~,grid_data.gen_bus_rotanti] = xlsread(filename, 'Generatori Rotanti_A', 'B2:B7');
grid_data.gen_P_rotanti = xlsread(filename, 'Generatori Rotanti_A', 'C2:C7');
grid_data.gen_Q_rotanti = xlsread(filename, 'Generatori Rotanti_A', 'D2:D7');
[~,grid_data.gen_type_rotanti] = xlsread(filename, 'Generatori Rotanti_A', 'H2:H7');
[~,grid_data.gen_rotanti_types_ins] = xlsread(filename, 'Generatori Rotanti_A', 'F2:F7');

[~,grid_data.gen_bus_statici] = xlsread(filename, 'Generatori Statici_A', 'B2:B23');
grid_data.gen_P_statici = xlsread(filename, 'Generatori Statici_A', 'C2:C23');
[~,grid_data.gen_statici_types_ins] = xlsread(filename, 'Generatori Statici_A', 'N2:N23');


[grid_data.gen_statici_nom,grid_data.gen_statici_types_all] = xlsread(filename2, 'Generatori Statici', 'A2:B20');
[grid_data.gen_rotanti_nom,grid_data.gen_rotanti_types_all] = xlsread(filename2, 'Generatori Rotanti', 'A2:B7');
