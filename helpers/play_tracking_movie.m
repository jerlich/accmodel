function M=play_tracking_movie(points,ts,x,y,theta,ptime)

if nargin<6
    ptime=0.03;
end

n_samp=size(points,2);

if nargin==1
    ts=1:n_samp;
end

[r,g,b]=vid_bitfield(points);

figure(747); clf
ax=axes('Position',[0.15 0.15 0.7 0.7]);
set(ax,'YLim',[0 700],'XLim',[0 700],'NextPlot','add');
for px=1:n_samp
    
    redidx=find(r(:,1)==px);
    redx=r(redidx,2);
    redy=r(redidx,3);

    blueidx=find(b(:,1)==px);
    bluex=b(blueidx,2);
    bluey=b(blueidx,3);

    greenidx=find(g(:,1)==px);
    greenx=g(greenidx,2);
    greeny=g(greenidx,3);

    plot(ax,redx,redy,'r.');
    plot(ax,greenx,greeny,'g.');
    plot(ax,bluex,bluey,'b.');
    text(20,20,sprintf('%5.3g',ts(px)));
    
    plot(ax,x(px),y(px),'+');
    plot(ax,x(px),y(px),'+');
    [u v]=pol2cart(theta(px)*pi/180,100);
    quiver(ax, x(px),y(px), u,v )
    
    M(px)=getframe;
    pause(ptime)
    cla;
    
end

    plot(ax,[10 500], [500 10],'k-');
    
    M(end+1)=getframe;

%movie(M);

%keyboard