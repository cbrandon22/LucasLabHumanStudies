function le_mm_testMotorSelection

subjList = {'HUP084','HUP087','HUP089','HUP090','HUP091','HUP092','HUP111'};
allMotor = [];
allNonmotor = [];
for i=1:length(subjList)
    subj=subjList{i}
    [eLbl_list_motor,eLbl_list_nonmotor] = le_mm_selectMotorSites (subj);
    allMotor = [allMotor,eLbl_list_motor];
    allNonmotor = [allNonmotor,eLbl_list_nonmotor];
end
keyboard;
