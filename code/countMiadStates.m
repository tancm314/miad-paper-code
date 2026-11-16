function [N_count,mean_dis_atti,mean_dis_vec]=countMiadStates(save_result,criterion)
%{
Revised by Caiming Tan on 16 December 2022 in Hongyunjiayuan community
Written by Caiming Tan on 8 December 2022 in Hongyunjiayuan community
%}
N=length(save_result);
N_wrong_selection=0;
N_no_alarm_at_wrong_selection=0;
N_right_selection=0;
N_alarm_at_right_selection=0;
save_dis_atti=[];
save_dis_vec=[];
for i=1:1:N
    dis_vec=save_result(i,1);
    dis_atti=save_result(i,2);
    alarm=save_result(i,3);
    if dis_vec>0
        save_dis_atti=[save_dis_atti,dis_atti*180/pi];
        save_dis_vec=[save_dis_vec,dis_vec];
        N_wrong_selection=N_wrong_selection+1;
        if alarm==false
            N_no_alarm_at_wrong_selection=N_no_alarm_at_wrong_selection+1;
        end
    else
        N_right_selection=N_right_selection+1;
        if alarm==true
            N_alarm_at_right_selection=N_alarm_at_right_selection+1;
        end
    end
end
N_count=[N_right_selection,N_alarm_at_right_selection,...
    N_wrong_selection,N_no_alarm_at_wrong_selection];

is_debug=false;
if N_wrong_selection
    mean_dis_atti=mean(save_dis_atti);
    mean_dis_vec=mean(save_dis_vec);
    if is_debug
        if criterion==1
            figure(2+criterion);
        else
            figure(2+criterion);
        end
        plot(save_dis_atti,'.-');
    end
else
    mean_dis_atti=0;
    mean_dis_vec=0;
end
