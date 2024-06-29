% OP data cope with
[OP_all]=xlsread('C:\Users\dell\Desktop\2023.0257\scr\output\OP_example_com_mean.xlsx');
OP_cope=OP_all./OP_all(:,1);
[OP_all_std]=xlsread('C:\Users\dell\Desktop\2023.0257\scr\output\OP_example_com_std.xlsx');
OP_cope_std=OP_all_std./OP_all_std(:,1);

xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\OP_example_compare_mean.xlsx',OP_cope);
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\OP_example_compare_std.xlsx',OP_cope_std);

% Wine data cope with

[Wine_all]=xlsread('C:\Users\dell\Desktop\2023.0257\scr\output\Wine_example_com_mean.xlsx');
Wine_cope=Wine_all./Wine_all(:,1);
[Wine_all_std]=xlsread('C:\Users\dell\Desktop\2023.0257\scr\output\Wine_example_com_std.xlsx');
Wine_cope_std=Wine_all_std./Wine_all_std(:,1);

xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Wine_example_compare_mean.xlsx',Wine_cope);
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Wine_example_compare_std.xlsx',Wine_cope_std);

% Piston data cope with

[Piston_all]=xlsread('C:\Users\dell\Desktop\2023.0257\scr\output\Piston_example_com_mean.xlsx');
Piston_cope=Piston_all./Piston_all(:,1);
[Piston_all_std]=xlsread('C:\Users\dell\Desktop\2023.0257\scr\output\Piston_example_com_std.xlsx');
Piston_cope_std=Piston_all_std./Piston_all_std(:,1);

xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Piston_example_compare_mean.xlsx',Piston_cope);
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\Piston_example_compare_std.xlsx',Piston_cope_std);
