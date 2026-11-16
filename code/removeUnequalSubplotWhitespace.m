% Revised by Caiming Tan on 1 November 2022 in Hongyunjiayuan community
 
% author : TSC
% time   : 2017-01-02
% email  : 292936085#qq.com(��#�滻Ϊ@)
function removeUnequalSubplotWhitespace(s_gca,sub_row,proportion)
N=length(s_gca);
s_pos=zeros(sub_row,4);
for i=1:1:sub_row
    gca_sub=s_gca(i);
    s_pos(i,:) = getAxesPositionWithoutWhitespace(gca_sub, sub_row, 1, i, 1);
end
pos_new_x=max(s_pos(:,1));
pos_new_w=min(s_pos(:,3));
pos_new_h_all=sum(s_pos(:,4));
s_pos(:,1)=pos_new_x;
s_pos(:,3)=pos_new_w;
s_pos_reverse=s_pos(end:-1:1,:);
s_dh=s_pos_reverse(2:end,2)-(s_pos_reverse(1:end-1,2)+s_pos_reverse(1:end-1,4));
for i=1:1:sub_row
    s_pos_reverse(i,4)=pos_new_h_all*proportion(sub_row-i+1)/sum(proportion);
end
for i=1:1:sub_row
    if i>1
        s_pos_reverse(i,2)=s_pos_reverse(i-1,2)+s_pos_reverse(i-1,4)+s_dh(i-1);
    end
end
s_pos=s_pos_reverse(end:-1:1,:);
for i=1:1:sub_row
    gca_sub=s_gca(i);
    set(gca_sub, 'Position', s_pos(i,:));
end
% ����Position
% set(gca, 'Position', [pos_new_x, pos_new_y, pos_new_w, pos_new_h]);

function pos = getAxesPositionWithoutWhitespace(gca, sub_row, sub_col, current_row, current_col)
% ����OuterPosition
sub_axes_x = current_col*1/sub_col - 1/sub_col;
sub_axes_y = 1-current_row*1/sub_row; % y�Ǵ������µ�
sub_axes_w = 1/sub_col;
sub_axes_h = 1/sub_row;
set(gca, 'OuterPosition', [sub_axes_x, sub_axes_y, sub_axes_w, sub_axes_h]); % ����OuterPosition
 
% TightInset��λ��
inset_vectior = get(gca, 'TightInset');
inset_x = inset_vectior(1);
inset_y = inset_vectior(2);
inset_w = inset_vectior(3);
inset_h = inset_vectior(4);
 
% OuterPosition��λ��
outer_vector = get(gca, 'OuterPosition');
pos_new_x = outer_vector(1) + inset_x; % ��Position��ԭ���Ƶ���TightInset��ԭ��
pos_new_y = outer_vector(2) + inset_y;
pos_new_w = outer_vector(3) - inset_w - inset_x; % ����Position�Ŀ�
pos_new_h = outer_vector(4) - inset_h - inset_y; % ����Position�ĸ�

pos=[pos_new_x, pos_new_y, pos_new_w, pos_new_h];