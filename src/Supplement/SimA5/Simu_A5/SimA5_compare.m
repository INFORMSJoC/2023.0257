[simA5_all]=xlsread('C:\Users\dell\Desktop\2023.0257\scr\output\SimA5_compare.xlsx');
simA5_com=simA5_all./simA5_all(:,3);
xlswrite('C:\Users\dell\Desktop\2023.0257\scr\output\SimA5_compare',simA5_com);