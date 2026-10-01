function Colors = f_Clrs_fresh(n, NW)
%F_CLRS_FRESH n distinct colors (RGB, n x 3).
%
%   Colors = f_Clrs_fresh(n)
%   Colors = f_Clrs_fresh(n, 'Order', Order)
%   Order: 'my' (default, sorted by hue), 'convenient' (fixed order),
%   'rainbow' (sorted by hue). For n > 20 the jet colormap is returned.
%   Colors are based on
%   https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
%
%   Alexey Ryabov

arguments
    n double
    NW.Order  {mustBeMember(NW.Order,{'my','convenient', 'rainbow'})} = 'my'
end

if n < 21
    switch NW.Order
        case 'my'
            Colors = [...
                204, 51, 0;
                0, 153, 255
                245, 130, 48;...% orange
                145, 30,  180;...% purple
                60,  180, 75;...% green
                170, 110, 40;...% brown
                255, 200, 25; % dark yellow
                230, 190, 255;...% lavender
                210, 245, 60;...% lime
                230, 25, 75;
                70, 240, 240;
                240, 50, 230;
                250, 190, 212;
                255, 250, 200;
                0,   128, 128;...% teal
                170, 255, 195;
                0, 0, 128;
                128, 128, 0;
                255, 215, 180;
                0,   130, 200;...% blue
                128, 128, 128;
                0, 0, 0
                ]/255;
            Colors = Colors(1:n, :);
            
            Colors = rgb2hsv(Colors);
            [~, ind] = sort(Colors(:, 1));
            Colors = Colors(ind, :);
            Colors = hsv2rgb(Colors);
        case 'convenient'
            Colors = [...
                230, 25, 75;
                60, 180, 75;...
                255, 200, 25; % dark yellow
                0, 0, 128; ...%Navi
                245, 130, 48;...% orange
                145, 30,  180;...% purple
                70, 240, 240;
                240, 50, 230; %Magenta
                210, 245, 60;...% lime
                250,190,212; ... %Pink
                128, 0, 0;... %Maroon
                230, 190, 255;...% lavender
                170, 110, 40;...% brown
                255, 250, 200;... %Beige
                0,   128, 128;...% teal
                170, 255, 195;...%Mint
                128, 128, 0;...%Olive
                255, 215, 180;... %Apricot
                0,   130, 200;...% blue
                128, 128, 128; ...%grey
                0, 0, 0
                ]/255;
            Colors = Colors(1:n, :);  
            case 'rainbow'
               Colors = [...
                230, 25, 75;
                60, 180, 75;...
                255, 200, 25; % dark yellow
                0,   130, 200;...% blue
                245, 130, 48;...% orange
                145, 30,  180;...% purple
                70, 240, 240;
                240, 50, 230; %Magenta
                210, 245, 60;...% lime
                250,190,212; ... %Pink
                0,   128, 128;...% teal
                230, 190, 255;...% lavender
                170, 110, 40;...% brown
                255, 250, 200;... %Beige
                128, 0, 0;... %Maroon
                170, 255, 195;...%Mint
                128, 128, 0;...%Olive
                255, 215, 180;... %Apricot
                0, 0, 128; ...%Navi
                ]/255;
            Colors = rgb2hsv(Colors);
            [~, ind] = sort(Colors(:, 3));
            Colors = Colors(ind(1:n), :);
             
            [~, ind] = sort(Colors(:, 1));
            Colors = Colors(ind, :);
            Colors = hsv2rgb(Colors);
    end
else
    Colors = jet(n);
end
