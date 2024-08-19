function time_out = time_addition_1day(time_in)

t = datetime(time_in.year,time_in.month,time_in.day);
t = t+day(1);
t2 = datevec(t);
time_out.year = t2(1);
time_out.month = t2(2);
time_out.day = t2(3);
time_out.week = weekday(t);

