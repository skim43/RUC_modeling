classdef T3DAssemble < T3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T3DAssemble class
% Using T3D class, assemble two T3D packs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 16-Feb-2020
% by Sangmin Lee, Ph.D.
% Research Scientist
% C-SWARM
% Center for Shock-Wave Processing of Advanced Reactive Materials
% 117 Cushing Hall, Notre Dame, IN 46556
% Department of Aerospace and Mechanical Engineering
% University of Notre Dame, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  properties
    list_1to2;
    list_2to1;
  end
  methods 
    function assemble_two_packs(obj, T1, T2)

      obj.E   = cell(obj.t3d.ids(end), 1);
      [obj.list_1to2, obj.list_2to1] = list_R_to_R(T1, T2, obj.tol);
      
      T1max = zeros(obj.t3d.ids(end), 1);
      for ia = 1: obj.t3d.ids(end)
        n1 = numel(T1.E{ia});
        n2 = numel(T2.E{ia});
        
        T1max(ia) = size(T1.map{ia}, 1);
        obj.E{ia}(1:n1)      = T1.E{ia}(:);
        obj.E{ia}(n1+(1:n2)) = T2.E{ia}(:);
      end
      
      % supress 2nd matrix
      obj.E{4}{end}.valid = false;
      
      % differnciate 2nd materials
      for ia = (numel(T1.E{4})+1): numel(obj.E{4})-1
        obj.E{4}{ia}.property = obj.E{4}{ia}.property + 1;
      end
      
      T1max = T1max + 100;
      % update inner particle ids
      for ia = 1: obj.t3d.ids(end)
        n1 = numel(T1.E{ia});        
        n2 = numel(T2.E{ia});
        switch(obj.t3d.entities{ia})
          case 'vertex'
            for ib = 1: n2
              obj.E{ia}{n1+ib}.id = obj.E{ia}{n1+ib}.id + T1max(ia);
            end
          case 'curve'
            for ib = 1: n2
              obj.E{ia}{n1+ib}.id = obj.E{ia}{n1+ib}.id + T1max(ia);
              obj.E{ia}{n1+ib}.e  = sign(obj.E{ia}{n1+ib}.e).*(abs(obj.E{ia}{n1+ib}.e) + T1max(1));
            end
          case 'surface'
            for ib = 1: n2
              obj.E{ia}{n1+ib}.id = obj.E{ia}{n1+ib}.id + T1max(ia);
              obj.E{ia}{n1+ib}.e  = sign(obj.E{ia}{n1+ib}.e).*(abs(obj.E{ia}{n1+ib}.e) + T1max(2));
            end
          case 'region'
            for ib = 1: n2
              obj.E{ia}{n1+ib}.id = obj.E{ia}{n1+ib}.id + T1max(ia);
              obj.E{ia}{n1+ib}.e{1}  = sign(obj.E{ia}{n1+ib}.e{1}).*(abs(obj.E{ia}{n1+ib}.e{1}) + T1max(3));
              obj.E{ia}{n1+ib}.e{2}  = sign(obj.E{ia}{n1+ib}.e{2}).*(abs(obj.E{ia}{n1+ib}.e{2}) + T1max(6));
              obj.E{ia}{n1+ib}.e{3}  = sign(obj.E{ia}{n1+ib}.e{3}).*(abs(obj.E{ia}{n1+ib}.e{3}) + T1max(5));
            end 
          case 'patch'
            for ib = 1: n2
              obj.E{ia}{n1+ib}.id = obj.E{ia}{n1+ib}.id + T1max(ia);
              obj.E{ia}{n1+ib}.e  = sign(obj.E{ia}{n1+ib}.e).*(abs(obj.E{ia}{n1+ib}.e) + T1max(2));
              if(numel(obj.E{ia}{n1+ib}.mirror)>0)
                obj.E{ia}{n1+ib}.mirror  = obj.E{ia}{n1+ib}.mirror + T1max(5);
              end
            end
          case 'shell'
            for ib = 1: n2
              obj.E{ia}{n1+ib}.id       = obj.E{ia}{n1+ib}.id + T1max(ia);
              obj.E{ia}{n1+ib}.e(1)     = obj.E{ia}{n1+ib}.e(1) +  T1max(3);
              obj.E{ia}{n1+ib}.e(2:end) = sign(obj.E{ia}{n1+ib}.e(2:end)).*(abs(obj.E{ia}{n1+ib}.e(2:end)) + T1max(2));
            end 
          case 'interface'
          
        end
      end
      
      % update map      
      obj.map = cell(obj.t3d.ids(end), 1);
 
      for ia = 1: obj.t3d.ids(end)
        limit = 0;
        for ib = 1: numel(obj.E{ia})
          limit = max(limit, obj.E{ia}{ib}.id);
        end
        obj.map{ia} = sparse(limit, 1);
      end
      
      for ia = 1: obj.t3d.ids(end)
        for ib = 1: numel(obj.E{ia})
          obj.map{ia}(obj.E{ia}{ib}.id) = ib;
        end
      end      

      obj.remove_duplicated_vertices;

      % minimum point distance
      % min_d = obj.compute_smallest_edge;
      
      obj.update_edges;
      obj.remove_patch_for_duplicated_matrx(T1, T2);
      obj.update_patch(T1, T2);
      obj.update_region(T1, T2);
    end
    function remove_duplicated_vertices(obj)
    %
      vno = numel(obj.E{1});
      vlist = zeros(vno, 1);
      for ia = 1: vno-1
        if(vlist(ia)>0); continue; end
        for ib = ia+1:vno
          if(vlist(ib)>0); continue; end
          if(norm(obj.E{1}{ia}.e - obj.E{1}{ib}.e)<obj.tol)
            vlist(ib) = ia;
            obj.E{1}{ib}.valid = false;
            obj.E{1}{ib}.id = obj.E{1}{ia}.id;
            if(numel(obj.E{1}{ia}.option)>2 && numel(obj.E{1}{ib}.option) < 2)
              obj.E{1}{ia}.option = obj.E{1}{ib}.option;
            end
          end
        end
      end
      for ia = 1: numel(obj.E{2})
        v = obj.E{2}{ia}.e;
        L = obj.map{1}(v);
        G = [obj.E{1}{L(1)}.id, obj.E{1}{L(2)}.id];
        obj.E{2}{ia}.e(vlist(L)>0) = G(vlist(L)>0);
      end
    end
    function plot_edges(obj)
      figure; hold on;

      for ia = 1: numel(obj.E{2})
        if(~obj.E{2}{ia}.valid || numel(obj.E{2}{ia}.order)>0); continue; end
        v = obj.map{1}(obj.E{2}{ia}.e);
        x = [obj.E{1}{v(1)}.e; obj.E{1}{v(2)}.e];
        plot3(x(:, 1), x(:, 2), x(:, 3), '-ok'); 
      end
      hold off;
      camlight
      xlabel('x');
      ylabel('y');
      zlabel('z');
      axis equal
      view(150,30);
      set(gcf,'Color','white');
    end
    function update_edges(obj)    
    % indentify out box
      xmax = zeros(1, 3) - 1.0e15;
      xmin = zeros(1, 3) + 1.0e15;
      v_on_ebnd = zeros(numel(obj.E{1}), 2);
      e_on_ebnd = zeros(numel(obj.E{2}), 2);

      for ia = 1: numel(obj.E{2})
        if(numel(obj.E{2}{ia}.order)>0); continue; end
        e_on_ebnd(ia, :) = [ia, obj.E{2}{ia}.id];
        v = obj.map{1}(obj.E{2}{ia}.e);
        for ib = 1: 2
          if(obj.E{1}{v(ib)}.valid)
            v_on_ebnd(v(ib), :) = [v(ib), obj.E{2}{ia}.e(ib)];
            xmax = max([xmax; obj.E{1}{v(ib, 1)}.e]);
            xmin = min([xmin; obj.E{1}{v(ib, 1)}.e]);
          end
        end
      end
      v_on_ebnd(v_on_ebnd(:, 1)==0, :) = [];
      e_on_ebnd(e_on_ebnd(:, 1)==0, :) = [];

      % construct out box edges
      box.v = [xmax(1), xmax(2), xmax(3)
               xmax(1), xmax(2), xmin(3)
               xmax(1), xmin(2), xmax(3)
               xmax(1), xmin(2), xmin(3)
               xmin(1), xmax(2), xmax(3)
               xmin(1), xmax(2), xmin(3)
               xmin(1), xmin(2), xmax(3)
               xmin(1), xmin(2), xmin(3)];

      box.e = [1 2; 3 4; 5 6; 7 8; 1 3; 2 4; 5 7; 6 8; 1 5; 2 6; 3 7; 4 8];

      box.eno = size(box.e, 1);
      box.n = zeros(box.eno, 3);
      for ia = 1: box.eno
        box.n( ia, :) = box.v(box.e(ia, 2), :) - box.v(box.e(ia, 1), :);
        box.mn(ia)    = norm(box.n(ia, :)); 
        box.n( ia, :) = box.n(ia, :)./box.mn(ia);
      end
            
      % collect vertices align on out box edges
      l = @(mm, p, qi) abs(norm(qi - p(1, :)) + norm(p(2, :) - qi) - mm) < obj.tol;
      vno = size(v_on_ebnd, 1);
      v = cell(box.eno, 1);
      new_eno = 0;
      for ia = 1: box.eno
        v{ia} = zeros(vno, 2);
        for ib = 1: vno
          if(l(box.mn(ia), box.v(box.e(ia, :), :), obj.E{1}{v_on_ebnd(ib, 1)}.e))
            v{ia}(ib, :) = [v_on_ebnd(ib, 2), norm(obj.E{1}{v_on_ebnd(ib, 1)}.e - box.v(box.e(ia, 1), :))];
          end
        end
        v{ia}(v{ia}(:, 1)==0, :) = [];
        [~, J] = sort(v{ia}(:, 2));
        v{ia} = v{ia}(J, :);
        new_eno = new_eno + size(v{ia}, 1) - 1;
      end

      % create new edges
      cno   = numel(obj.E{2});
      mapno = numel(obj.map{2});
      obj.map{2}(end + (1:new_eno)) = zeros(new_eno, 1);
      obj.E{  2}(end + (1:new_eno)) =  cell(new_eno, 1);
      
      cnt = 0;
      c0  = obj.E{2}{e_on_ebnd(1, 1)}; 
      for ia = 1: box.eno
        for ib = 1: (size(v{ia}, 1) - 1)
          cnt = cnt + 1;
          obj.E{2}{cno + cnt} = c0;
          obj.E{2}{cno + cnt}.id = mapno + cnt;
          obj.E{2}{cno + cnt}.e  = [v{ia}(ib, 1), v{ia}(ib+1, 1)];
          %obj.E{1}{obj.map{1}(v{ia}(ib,   1))}.option = ''; % supress virtual vertex
          %obj.E{1}{obj.map{1}(v{ia}(ib+1, 1))}.option = ''; % supress virtual vertex         
          obj.map{2}(mapno + cnt) = cno + cnt;
        end
      end


      % disable edges for replacing with new edges
      for ia = 1: size(e_on_ebnd, 1)
        obj.E{2}{e_on_ebnd(ia, 1)}.valid = false;
      end
      
      % create list of replace old edges to new edges
      old2new = sparse(mapno+new_eno, 1);
      replace = cell(size(e_on_ebnd, 1), 2);
      for ia = 1: size(e_on_ebnd, 1)
        tmp = obj.E{2}{e_on_ebnd(ia, 1)}.e;
        for ib = 1: box.eno
          ck  = zeros(1,2);
          for ic = 1: size(v{ib}, 1)
            if(tmp(1)==v{ib}(ic, 1))
              ck(1) = ic;
            end
            if(tmp(2)==v{ib}(ic, 1))
              ck(2) = ic;
            end
          end
          if(ck(1)>0 && ck(2)>0)
            if(ck(2)>ck(1))
              m = v{ib}(ck(1):ck(2), 1);
            else
              m = v{ib}(ck(1):-1:ck(2), 1);
            end
            replace{ia, 1} = [m(1:end-1), m(2:end)];
            old2new(e_on_ebnd(ia, 2)) = ia;
            break;
          end          
        end
      end
      
      % assign new edge ids to the replacements
      for ia = 1: size(replace, 1)
        replace{ia, 2} = zeros(1, size(replace{ia, 1}, 1));
        for ib = 1: size(replace{ia, 1}, 1)
          tmp = replace{ia, 1}(ib, :);
          for ic = (cno+1):1:(cno+new_eno)
            if(tmp(1) == obj.E{2}{ic}.e(1) && tmp(2) == obj.E{2}{ic}.e(2))
              replace{ia, 2}(ib) = obj.E{2}{ic}.id;
              break;
            end

            if(tmp(2) == obj.E{2}{ic}.e(1) && tmp(1) == obj.E{2}{ic}.e(2))
              replace{ia, 2}(ib) = -obj.E{2}{ic}.id;
              break;
            end            
          end
        end
      end      
      
      % update face and shell and surface
      surf = [3, 5, 6];
      for iA = 1: numel(surf)
        for ia = 1: numel(obj.E{surf(iA)})
          if(surf(iA)==6)
            c = obj.E{surf(iA)}{ia}.e(2:end);
          else
            c = obj.E{surf(iA)}{ia}.e;
          end
          is_updated = false;
          for ib = numel(c):-1:1
            o2n = old2new(abs(c(ib)));
            if(o2n>0)
              is_updated = true;

              if(c(ib)<0)
                c = [c(1:(ib-1)), -replace{o2n, 2}(end:-1:1), c((ib+1): end)];
              else
                c = [c(1:(ib-1)),  replace{o2n, 2},           c((ib+1): end)];
              end 
            end
          end
          if(is_updated)
            if(surf(iA)==6)
              obj.E{surf(iA)}{ia}.e = [obj.E{surf(iA)}{ia}.e(1), c]; % 1st entity of shell is background surface id
            else
              obj.E{surf(iA)}{ia}.e = c;
            end
          end
        end
      end
    end
    function remove_patch_for_duplicated_matrx(obj, T1, T2)    
      % disable patches used only for the matrix for T2
      pno = numel(T1.E{5});
      for ia = 1: numel(T2.E{5})
        obj.E{5}{pno + ia}.valid = false;        
      end
      
      rno = numel(T1.E{4});
      for ia = 1: numel(T2.E{4})-1
        for ib = 1: numel(obj.E{4}{rno + ia}.e{3})
          obj.E{5}{obj.map{5}(abs(obj.E{4}{rno + ia}.e{3}(ib)))}.valid = true;
        end        
      end
    end
    function update_patch(obj, T1, T2)
      % substract T2 patches from T1 on out box
      rno = numel(T1.E{4});
      for iA = 1: (rno-1)
        iB = obj.list_1to2(iA);
        if(iB==0); continue; end
        f1 = obj.E{4}{iA    }.e{3};
        f2 = obj.E{4}{iB+rno}.e{3};
        if(numel(f1)==0 || numel(f2)==0); continue; end        
        for ia = 1: numel(f1)
          m = obj.E{5}{obj.map{5}(abs(f1(ia)))}.n;
          m = m./norm(m);
          c_ia = obj.E{5}{obj.map{5}(abs(f1(ia)))}.e;
          for ib = 1: numel(f2)
            n = obj.E{5}{obj.map{5}(abs(f2(ib)))}.n;
            n = n./norm(n);            
            c_ib = obj.E{5}{obj.map{5}(abs(f2(ib)))}.e;
            
            if(abs(m*n(:) - 1)<obj.tol)
              tmp = [abs(c_ia), abs(c_ib); 1:numel(c_ia), 1:numel(c_ib)]';
              [~, ix, jx] = unique(tmp(:, 1));

              vec=histc(jx,1:max(jx));
              qx = vec==1;
              tmp(ix(qx), :) = [];
              if(size(tmp, 1)>0)
                if(size(tmp, 1) == 4)
                  c_ia(tmp(2,2)) = [];
                  c_ib(tmp(4,2)) = [];
                  tmp([2,4], :)  = [];
                end
                if(sign(c_ia(tmp(1,2))) ~= sign(c_ib(tmp(2,2))))
                  obj.E{5}{obj.map{5}(abs(f1(ia)))}.e = [c_ia(1:tmp(1,2)-1), ...
                                                         c_ib(tmp(2,2)+1:end), ...
                                                         c_ib(1:(tmp(2,2)-1)), ...
                                                         c_ia(tmp(1,2)+1: end)];                  
                else
                  obj.E{5}{obj.map{5}(abs(f1(ia)))}.e = [c_ia(1:tmp(1,2)-1), ...
                                                         -c_ib((tmp(2,2)-1):-1:1), ...
                                                         -c_ib(end:-1:(tmp(2,2)+1)), ...
                                                         c_ia(tmp(1,2)+1: end)];
                end
              else
                if(sign(m*n(:))<0)
                  obj.E{5}{obj.map{5}(abs(f1(ia)))}.e = [c_ia, c_ib];
                else
                  obj.E{5}{obj.map{5}(abs(f1(ia)))}.e = [c_ia, -c_ib(end:-1:1)];
                end
              end
              break;               
            end
          end
        end
      end      
    end
    
    function update_region(obj, T1, T2)
      % conform regions
      rno = numel(T1.E{4});
      surf = [1, 2];
      for iA = 1: (rno-1)
        iB = obj.list_1to2(iA);
        if(iB == 0); continue; end
        for ia = 1: numel(surf)
          f1 = obj.E{4}{iA    }.e{surf(ia)};
          f2 = obj.E{4}{iB+rno}.e{surf(ia)};
          if(numel(f2)==0); continue; end
          obj.E{4}{iA}.e{surf(ia)} = [f1, -f2];
        end
      end
      % remove vertual option for physical boundary
      surf = [3,5];
      id   = [1,3];
      for iA = 1: numel(surf)
        f = obj.E{4}{rno}.e{id(iA)};
        for ia = 1: numel(f)
          obj.E{surf(iA)}{obj.map{surf(iA)}(abs(f(ia)))}.option = '';
          c = obj.E{surf(iA)}{obj.map{surf(iA)}(abs(f(ia)))}.e;
          for ib = 1: numel(c)
            obj.E{2}{obj.map{2}(abs(c(ib)))}.option = '';
            v = obj.E{2}{obj.map{2}(abs(c(ib)))}.e;
            for ic = 1: numel(v)
              obj.E{1}{obj.map{1}(abs(v(ic)))}.option = '';
            end
          end
        end
      end      
    end
  end
end

function [list_1to2, list_2to1] = list_R_to_R(T1, T2, tol)
  rno1 = numel(T1.E{4});
  rno2 = numel(T2.E{4});
  list_1to2 = zeros(rno1, 1);
  list_2to1 = zeros(rno2, 1);

  for ia = 1: (rno1 - 1)
    cp1 = compute_cp(T1, ia);
    if(numel(cp1)==0); continue; end
    for ib = 1: (rno2 - 1)
      if(list_2to1(ib)>0); continue; end;
      cp2 = compute_cp(T2, ib);
      if(numel(cp1)==0); continue; end
      if(norm(cp1 - cp2)<tol)
        list_1to2(ia) = ib;
        list_2to1(ib) = ia;
        break;
      end
    end
  end
end

function cp = compute_cp(T, id)
  % center point of a sphere (T.E{4}{id})
  cp = [];

  surf = [3, 6];
  f = T.E{4}{id}.e;
  if(numel(f{1}) == 0 && numel(f{2})==0)
    return;
  end
  for ia = 1: 2
    if(numel(f{ia})==0); continue; end
    if(surf(ia)==6) % if shell
      sh = f{ia}(1);
      s  = T.E{6}{T.map{6}(abs(sh))}.e(1);
    else
      s = f{ia}(1);
    end
    c = T.E{3}{T.map{3}(abs(s))}.e(1);
    v = T.E{2}{T.map{2}(abs(c))}.e;
    x = [T.E{1}{T.map{1}(abs(v(1)))}.e; T.E{1}{T.map{1}(abs(v(2)))}.e];
    
    cp = mean(x);
    break;
  end
end