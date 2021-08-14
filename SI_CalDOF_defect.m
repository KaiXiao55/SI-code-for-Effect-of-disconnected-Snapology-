%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mobility of 3D tesselation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all;
% ------------------------------------------------- Input   ------------------------------------------------- 
tessellateType='Square';  
defectCase= 'a.2';    % a.1 represents the case in Figure 8 of the main text 
ijk=[2 2 2];  % tessellation pattern

 % ---------------------------------------------------------------------------------------------------------
Face=[3 4 2 1;1 2 6 5;2 4 8 6;5 7 3 1;7 8 4 3;5 6 8 7];
Edge=[1 2;1 3;3 4;2 4;1 5;2 6;3 7;4 8;5 6;6 8;5 7;7 8];
a1a=90*pi/180;  a2a=90*pi/180; a3a=90*pi/180;
beta=acos((cos(a1a)*cos(a2a)-cos(a3a))/(sin(a1a)*sin(a2a)))*180/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Main function %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pinitial=[0,0,0];

if strcmp(tessellateType,'Square')==1
[Node,refNode,PositionVec]=chaCoorWithRef(a1a,a2a,a3a,Pinitial);    %%unit cell   %%%%chaCoorWithRef(r1,r2,r3,r4,r5,Pinitial)要改
[extrudedFace,Node,extrudedEdge,exHinge]=extrudeNode(Node,Face,Edge);  %extrude unit cell
locLaV(1,:)=Node(19,:)-Node(22,:);    locLaV(2,:)=Node(26,:)-Node(15,:);    locLaV(3,:)=Node(31,:)-Node(10,:);

[mnode,mface,medge,mhinge,viewpoint]=ConstructDefectUnit(Node,extrudedFace,extrudedEdge,exHinge,locLaV,defectCase,[]);
elseif  strcmp(tessellateType,'Hex')==1
     [mnode,medge,mface,locLaV,viewpoint]=ConstructDefectHex('full');
end

figure(1)
PlotShape(mface,mnode,1,'txext',viewpoint)
mobilityUnit=CalMobility(mnode,mface,medge);
str1=('Mobility = ');
str2= num2str(mobilityUnit);
SC=[str1,str2];
annotation('textbox', [0.1, 0, 0.8, 0.1], 'String', SC,'fontsize',24,'linestyle','none','horizontalalignment','center');

%%-------------------------------
% Tessellate
gloLaV=2*locLaV;
[Asmnode,Asmface,Asmedge,Ashinge,viewpoint0]=ConstructDefectUnit(mnode,mface,medge,mhinge,gloLaV,'asssmbly',ijk);
figure(2)
PlotShape(Asmface,Asmnode,1,'tnxt',viewpoint)   % original setting
mobility=CalMobility(Asmnode,Asmface,Asmedge);
str1=('Mobility = ');
str2= num2str(mobility);
SC=[str1,str2];
annotation('textbox', [0.1, 0, 0.8, 0.1], 'String', SC,'fontsize',24,'linestyle','none','horizontalalignment','center');


%%%%% sub-function %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mnode,medge,mface,locLaV,viewpoint]=ConstructDefectHex(defectCase)
refNode=[0 0 0;  -1,-1.73205080756888,0;   -2,0,0;    0,0,3;];
switch defectCase
    case 'full'
        Pinitial=[refNode(1,:); refNode(2,:); refNode(3,:);   refNode(1,:)+  refNode(4,:);   refNode(2,:)+  refNode(4,:) ;   refNode(3,:)+ refNode(4,:);  ];
        index=[1 2 3 1 2 3 ];
        viewpoint=[145,34];
    case 'Planefull'
        Pinitial=[refNode(1,:); refNode(2,:); refNode(3,:);    ];
        index=[1 2 3  ];
        viewpoint=[145,34];
    case 'c0n4'
        Pinitial=[refNode(1,:); refNode(2,:);  refNode(2,:)+  refNode(4,:) ;   refNode(3,:)+ refNode(4,:);  ];
        index=[1 2  2 3  ];
        viewpoint=[145,34];
    case 'c1n4'
        Pinitial=[refNode(1,:);  refNode(2,:);  refNode(3,:);   ...
             refNode(2,:)+  refNode(4,:) ;  ];
        index=[1 2  3  2];
        viewpoint=[145,34];
    case 'c2n5'
        Pinitial=[refNode(1,:); refNode(2,:); refNode(3,:);    refNode(2,:)+  refNode(4,:) ;   refNode(3,:)+ refNode(4,:);  ];
          index=[1 2 3  2 3  ];
        viewpoint=[145,34];
end
unitCellFace=[];   
unitCellNode=[];  nNode=32; 
unitCellEdge=[];
for i=1:length(index)
    [Face{i},Node{i},Edge{i}]=ConstructHexPoly(90*pi/180,90*pi/180,90*pi/180,60*pi/180,Pinitial(i,:),index(i));
    unitCellNode(nNode*(i-1)+1:nNode*i,:)=Node{i};
    unitCellFace=[unitCellFace;  Face{i}+(i-1)*nNode];
    unitCellEdge=[unitCellEdge;  Edge{i}+(i-1)*nNode];
end
[mnode,medge,mface]=mergeNodeAdvance(unitCellNode,unitCellEdge,unitCellFace);
locLaV=[4+sqrt(3) 0 0 ; 2+sqrt(3)/2  -(2+sqrt(3)/2)*sqrt(3) 0; 0 0 6; ]/2;
end

function [extrudedFace,extrudedNode,extrudedEdge]=ConstructHexPoly(a1,a2,a3,b,Pinitial,index)  
L=1; 
V11=[0,0,L]; x1=[L,0,0]; 
z1=[0,0,L];  y1dot=[L*cos(b),L*sin(b),0];  x1dot=cross(y1dot,z1);
V12=cos(a1-pi/2)*x1dot+cos(a1)*z1;  V13=[0,L*cos(a2-pi/2),cos(a2)];

V21=z1;
x2=-y1dot;  z2=z1;  y2=cross(z2,x2);
V23=L*cos(a1-pi/2)*y2+L*cos(a1)*z2;
y2dot=-y1dot+x1; x2dot=cross(y2dot,z2); z2dot=z1;
V22=cos(a3-pi/2)*x2dot+cos(a3)*z2dot;

V31=z1;
x3=-y2dot; z3=z1; y3=cross(z3,x3);
V33=L*cos(a3-pi/2)*y3+L*cos(a3)*z3;
y3dot=-x1; x3dot=cross(y3dot,z3); z3dot=z3;
V32=cos(a2-pi/2)*x3dot+cos(a2)*z3dot;

switch index
    case 1
        P1=V12; P2=V13; P3=V11;
    case 2
        P1=V22; P2=V23; P3=V21;
    case 3
        P1=V32; P2=V33; P3=V31;
end
node1=Pinitial; node2=node1+P1; node3=node1+P2; node4=node1+P1+P2; node5=node1+P3; 
node6=node1+P1+P3; node7=node3+P3; node8=node4+P3;
ArchiNode=[real(node1);real(node2);real(node3);real(node4);real(node5);real(node6);real(node7);real(node8)];
Face=[3 4 2 1;1 2 6 5;2 4 8 6;5 7 3 1;7 8 4 3;5 6 8 7];
Edge=[1 2;1 3;3 4;2 4;1 5;2 6;3 7;4 8;5 6;6 8;5 7;7 8];
[extrudedFace,extrudedNode,extrudedEdge]=extrudeNode(ArchiNode,Face,Edge);
end

function [mnode,mface,medge,mhinge,viewpoint]=ConstructDefectUnit(unitNode,unitFace,unitEdge,unitHinge,lav,defectCase,layout)
viewpoint=[];
switch defectCase
    case 'planeFull'
        layPos=[[0 0 0]; lav(1,:);lav(2,:); lav(1,:)+lav(2,:); ];
        viewpoint=[127,44];
    case 'planePartial'
        layPos=[[0 0 0]; lav(1,:);lav(2,:)+lav(1,:); ];
        viewpoint=[127,44];
   case 'a.1'
        layPos=[[0 0 0]; lav(1,:);lav(2,:);lav(1,:)+lav(2,:);  lav(1,:)+lav(3,:); lav(3,:); lav(1,:)+lav(2,:)+lav(3,:); lav(2,:)+lav(3,:);];
        viewpoint=[142,20];
    case 'a.2'
        layPos=[[0 0 0]; lav(1,:);lav(2,:);lav(3,:); ];
        viewpoint=[142,20];
    case 'a.3'
        layPos=[[0 0 0]; lav(1,:);lav(2,:); lav(1,:)+lav(3,:); ];
        viewpoint=[142,20];
    case 'a.4'
        layPos=[[0 0 0]; lav(1,:);lav(2,:); lav(1,:)+lav(3,:); lav(2,:)+lav(3,:);];
        viewpoint=[142,20];
    case 'a.5'
        layPos=[[0 0 0]; lav(1,:);lav(2,:); lav(1,:)+lav(2,:);lav(3,:);];
        viewpoint=[142,20];
%           viewpoint=[119,17];
    case 'a.6'
        layPos=[[0 0 0]; lav(1,:);lav(2,:);lav(1,:)+lav(2,:);  lav(1,:)+lav(3,:); lav(2,:)+lav(3,:);];
        viewpoint=[142,20];
%           viewpoint=[125,16];
    case 'a.7'
        layPos=[[0 0 0]; lav(1,:);lav(2,:);lav(1,:)+lav(2,:);  lav(1,:)+lav(3,:); lav(3,:);];
        viewpoint=[142,20];
%          viewpoint=[152,13];
    case 'a.8'
        layPos=[[0 0 0]; lav(1,:);lav(2,:);lav(1,:)+lav(2,:);  lav(1,:)+lav(3,:); lav(3,:); lav(3,:)+lav(2,:);];
        viewpoint=[142,20];
%         viewpoint=[-130,20];
    case 'a.9'
        layPos=[ lav(1,:);lav(2,:);lav(1,:)+lav(2,:);  lav(1,:)+lav(3,:); lav(3,:); lav(3,:)+lav(2,:);];
        viewpoint=[142,20];
%         viewpoint=[140,16];
    case 'asssmbly'
        [layPos]=findLayoutPosition(lav,layout);
end
% Reconstruct the node
nNode=size(unitNode,1); face=[]; edge=[];  hinge=[];
for i=1:size(layPos,1)
    node(nNode*(i-1)+1:nNode*i,:)=unitNode+layPos(i,:);
    face=[face;unitFace+(i-1)*nNode];
    edge=[edge;unitEdge+(i-1)*nNode];
    hinge=[hinge;unitHinge+(i-1)*nNode];
end
[mnode,medge,mface,mhinge]=mergeNodeAdvance(node,edge,face,hinge);
end

function [layPos]=findLayoutPosition(Lav,layout)
m=1;
for k=1:layout(3)
    for j=1:layout(2)
        for i=1:layout(1)
            layPos(m,:)=[0 0 0]+(i-1)*Lav(1,:)+(j-1)*Lav(2,:)+(k-1)*Lav(3,:);
            m=m+1;
        end
    end
end
end

function [Node,refNode,PositionVec]=chaCoorWithRef(a1a,a2a,a3a,Pinitial)
L=1; e1=[1,0,0]; e2=[0,1,0]; e3=[0,0,1];
beta=acos((cos(a1a)*cos(a2a)-cos(a3a))/(sin(a1a)*sin(a2a)));
P1=L*sin(a2a)*sin(beta)*e1-L*sin(a2a)*cos(beta)*e2+L*cos(a2a)*e3;
P2=L*sin(a1a)*e2+L*cos(a1a)*e3;
P3=L*e3;

node1=Pinitial; node2=node1+P1; node3=node1+P2; node4=node1+P1+P2; node5=node1+P3;
node6=node1+P1+P3; node7=node3+P3; node8=node4+P3;
Node=[real(node1);real(node2);real(node3);real(node4);real(node5);real(node6);real(node7);real(node8)];

PositionVec=[];
PositionVec(1,:)=2*cross(P2,P3)/sqrt((P2*P2')*(P3*P3')-(P2*P3')^2);
PositionVec(2,:)=2*cross(P3,P1)/sqrt((P3*P3')*(P1*P1')-(P3*P1')^2);
PositionVec(3,:)=2*cross(P1,P2)/sqrt((P1*P1')*(P2*P2')-(P1*P2')^2);

refNode=[];
refNode(1,:)=real(node2)+PositionVec(1,:);  %%  (not use)
refNode(2,:)=real(node1)-PositionVec(2,:);  %%% second
refNode(3,:)=real(node1)-PositionVec(1,:)-PositionVec(2,:); %%% third
refNode(4,:)=real(node1)-PositionVec(1,:);   %%% Fourth
end

%%%%change number of node for other polyhedron



%%% Develop a function for extruding unit cell
function [extrudedFace,Node,extrudedEdge,exHinge]=extrudeNode(Node,Face,Edge)
k=size(Node,1)+1; e=size(Edge,1)+1; extrudedFace=zeros(4*size(Face,1),size(Face,2)); extrudedEdge=Edge; m=1;
exHinge=[];
for i=1:size(Face,1)
    basVec1=Node(Face(i,2),:)-Node(Face(i,1),:);  basVec2=Node(Face(i,3),:)-Node(Face(i,1),:);  %%%保证face中的点的编号都是顺时针，这样计算法向的refvector的方向才会正确
    refVec=cross(basVec1,basVec2)/sqrt((basVec1*basVec1')*(basVec2*basVec2')-(basVec1*basVec2')^2); %%the reference vector normal to the face for extruding cube
    for j=1:size(Face,2)
        Node(k,:)=real(Node(Face(i,j),:)+refVec); %%finished the extruding process of node
        if j<size(Face,2)
            extrudedFace(m,:)=[Face(i,j),Face(i,j+1),k+1,k];
            extrudedEdge(e,:)=[Face(i,j),k];
            extrudedEdge(e+1,:)=[k,k+1];
            extrudedEdge(e+2,:)=[Face(i,j),k+1];
            if j==1
            exHinge(end+1,:)=[ k  Face(i,1)  Face(i,j+1) Face(i,end)];
            elseif j<=size(Face,2)-1
                exHinge(end+1,:)=[ k   Face(i,j)  Face(i,j+1) Face(i,j-1)];
            end
            k=k+1;
            m=m+1;
            e=e+3;
        else
            extrudedFace(m,:)=[Face(i,j+1-size(Face,2)),Face(i,j),k,k+1-size(Face,2)]; %%the situation that j is equal to the size(Face,2)
            extrudedEdge(e,:)=[Face(i,j),k];
            extrudedEdge(e+1,:)=[k,k+1-size(Face,2)];
            extrudedEdge(e+2,:)=[Face(i,j),k+1-size(Face,2)];
            exHinge(end+1,:)=[ k   Face(i,j)   Face(i,1) Face(i,j-1)];
            k=k+1;
            m=m+1;
            e=e+3;
        end
    end
end
    for jex=1:size(Edge,1)
        kex=1;
        for iex=1:size(extrudedFace,1)
            if ismember(Edge(jex,1), extrudedFace(iex,:))==1 && ismember(Edge(jex,2), extrudedFace(iex,:))==1
                candiFace=extrudedFace(iex,:);
                candiFace(candiFace==Edge(jex,1))=[];
                candiFace(candiFace==Edge(jex,2))=[];
                hingeNode(kex)=candiFace(1);
                kex=kex+1;
            end
        end
         exHinge(end+1,:)=[Edge(jex,1) Edge(jex, 2) hingeNode(1) hingeNode(2)];
    end
end

function PlotShape(a,b,i,Text,viewpoint)
%clf;
extrudedFace=a;
Node=b;
patch('Faces',extrudedFace,'Vertices',Node,'FaceColor','m','facealpha',0.5); %this method is not appropriate for those polyhedrons whose faces have uncertain vertices
hold on ;
if strcmp(Text,'text')==1
    str={};
    for i=1:10000
        s=num2str(i);
        str{i}=s;
    end
    for j=1:size(Node,1)
        text(Node(j,1),Node(j,2),Node(j,3),str{j},'fontsize',10);
    end
end
view (viewpoint(1),viewpoint(2))%(125,37) % Feb. 16 %(120,59)%(149,14)% Nov.2rd
axis equal
axis off
light('color','w','style','infinite','position',3*[-0.5,1,0.5])
lighting flat
material dull
set(gcf,'color','w')
xlabel('x'),ylabel('y'),zlabel('z')  % 12.13 it's strange for cube that it doesn't transform the shape after transforming.
end

function [mnode,medge,mface,mhinge]=mergeNodeAdvance(node,edge,face,hinge)
% This version will delete the overplapping nodes and update the
% information in the edge and face
nM=1; % the cell matrix index colloect the information of overlapping node
for a=1:size(node,1)
    for b=1:size(node,1)
        if a~=b
            distance=sqrt((node(b,1)-node(a,1))^2+(node(b,2)-node(a,2))^2+(node(b,3)-node(a,3))^2);
            if distance<=10e-4
                mNode(nM,:)=[a b];
                nM=nM+1;
            end
        end
    end
end

mNodeUp=unique(sort(mNode,2),'rows');  % Keep the key information of deleteing nodes, the second column is the nodes needed to be deleted
% Keep only one nodes when one position has more than two nodes
deRo=[];
for i=1:size(mNodeUp,1)
    if ismember(mNodeUp(i,1),mNodeUp(:,2))==1
        if ismember(mNodeUp(i,2),mNodeUp(:,2))==1
            deRo=[deRo;i;];
        end
    end
end
mNodeUp(deRo,:)=[];
slave=mNodeUp(:,2);
% Replace the index of slave overlapped nodes
edgeUp=edge;
faceUp=face;
hingeUp=hinge;
for i=1:size(slave,1)
    edgeUp(edge==slave(i))=mNodeUp(i,1);  % update the element in the edges
    for j=1:size(face,1)
        faceUp(j,face(j,:)==slave(i))=mNodeUp(i,1);
    end
    for k=1:size(hinge,1)
        hingeUp(k,hinge(k,:)==slave(i))=mNodeUp(i,1);
    end
end
edgeUp=unique(sort(edgeUp,2),'rows');  % delete redundant edge constraints, but did not find a way to effiently deleted the

% Delete overlapping node
nindex = linspace(1,size(node,1),size(node,1));
nindex(mNodeUp(:,2))=[];
node(mNodeUp(:,2),:)=[];

% Update the new index to the edge and face
edgeUp2=edgeUp;
faceUp2=faceUp;
hingeUp2=hingeUp;
for i=1:size(nindex,2)
    edgeUp2(edgeUp==nindex(i))=i;  % update the element in the edges
    for j=1:size(faceUp,1)
        faceUp2(j,faceUp(j,:)==nindex(i))=i;
    end
    for k=1:size(hingeUp,1)
        hingeUp2(k, hingeUp(k,:)==nindex(i))=i;
    end
end
edgeUp2=unique(sort(edgeUp2,2),'rows');
mnode=node;  medge=edgeUp2;  mface=faceUp2;    mhinge=hingeUp2;% mface is not uniqued
end

%%%% Calculation of the mobility of entire structure
function mobility=CalMobility(Node,extrudedFace1,extrudedEdge1)
node=Node;
face=extrudedFace1;
edge=extrudedEdge1;%[1,2;2,3;3,4;4,1;1,5;2,6;3,7;3,8;4,9;1,10;5,11;6,12;7,13;8,13;9,14;10,11;11,12;12,13;13,14;14,11;1,15;4,16;3,17;2,18;18,15;18,1;15,16;15,4;16,17;16,3;17,18;17,2;1,19;2,20;6,21;5,22;22,19;22,1;19,20;19,2;20,21;20,6;21,22;21,5;3,23;7,24;6,25;2,26;26,23;26,3;23,24;23,7;24,25;24,6;25,26;25,2;3,27;8,28;13,29;7,30;30,27;30,3;27,28;27,8;28,29;28,13;29,30;29,7;3,31;4,32;9,33;8,34;34,31;34,3;31,32;31,4;32,33;32,9;33,34;33,8;4,35;1,36;10,37;9,38;38,35;38,4;35,36;35,1;36,37;36,10;37,38;37,9;1,39;5,40;11,41;10,42;42,39;42,1;39,40;39,5;40,41;40,11;41,42;41,10;5,43;6,44;12,45;11,46;46,43;46,5;43,44;43,6;44,45;44,12;45,46;45,11;6,47;7,48;13,49;12,50;50,47;50,6;47,48;47,7;48,49;48,13;49,50;49,12;8,51;9,52;14,53;13,54;54,51;54,8;51,52;51,9;52,53;52,14;53,54;53,13;9,55;10,56;11,57;14,58;58,55;58,9;55,56;55,10;56,57;56,11;57,58;57,14;11,59;12,60;13,61;14,62;62,59;62,11;59,60;59,12;60,61;60,13;61,62;61,14];
%%extrudedEdge1;
%%%%calcalate mobility
numEdge=size(edge,1); %number of the edges of the polyhedra
numFace=size(face,1); %number of the faces of the polyhedra
numNode=size(node,1)*3; %number of the nodes of the polyhedra*3
ja=[];
jw=[];
%length constraints
for j=1:numEdge
    for i=1:3
        k(j,i)=node(edge(j,1),i)-node(edge(j,2),i);
        jw(j,(edge(j,1)-1)*3+i)=k(j,i);
        jw(j,(edge(j,2)-1)*3+i)=-k(j,i);
    end
end
% n=size(jw,2)-rank(jw)-6; %  DOF of PFWs

%平面约束条件，对于cube计算不出ja正确的秩，但是对于hexagonal prism能计算出正确的秩
num=size(face,1);%确定平面的数量，以及循环的次数
x=1; %x为Ja矩阵的行数

for j=1:num      %单个平面的循环
    if (size(face(j,:),2)>3)%组成平面的点大于等于4个的情况，才能构成约束矩阵JA
        for r=4:size(face(j,:),2) %r为平面基点123之后的点
            pf1=node(face(j,1),:);
            pf2=node(face(j,2),:);
            pf3=node(face(j,3),:);
            pf4=node(face(j,r),:);
            kf1=cross((pf3-pf2),(pf3-pf4));
            kf2=cross((pf3-pf1),(pf4-pf1));
            kf3=cross((pf4-pf1),(pf2-pf1));
            kf4=cross((pf2-pf1),(pf3-pf1));
            for i=1:3
                ja(x,(face(j,1)-1)*3+i)=kf1(i);
                ja(x,(face(j,2)-1)*3+i)=kf2(i);
                ja(x,(face(j,3)-1)*3+i)=kf3(i);
                ja(x,(face(j,r)-1)*3+i)=kf4(i);
            end
            x=x+1;  %行数加一,已经添加了一个约束条件
        end
    end  %判断语句结束
end  %已经建立了一个平面的约束条件
%Constrains globe matrice
jj=[];
for i=1:numEdge
    for j=1:numNode
        jj(i,j)=jw(i,j);
    end
end
for i=1:x-1
    for j=1:numNode
        jj(i+numEdge,j)=ja(i,j);
    end
end
mobility=size(jj,2)-rank(jj)-6; % this equation is correct
% mobility=size(null(jj),2)-6;  %this equation is also correct
% ji2=pinv(jj); %jj is the jacobian matrix
% jt=(1-ji2*jj);
% DOF1=size(null(jw),2)-6;
end


