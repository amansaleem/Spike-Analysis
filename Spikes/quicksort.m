function[inData] = quicksort(inData, i, j, member)
% function[inData] = quicksort(inData, i, j, member)
% applies quicksort sorting algorithm to struct
% i is the starting index
% j i sthe end index
% member is the struct.member according to which struct should be sorted

% Oct-30 2007 LB

[pivot, pivotFound] =  findPivot(inData, i, j, member);
if pivotFound == 1
    [inData, p] = partition(inData, i, j, pivot, member);
    inData = quicksort(inData, i, p - 1, member);
    inData = quicksort(inData, p, j, member);
end

%--------------------------------------
function[pivot, pivotFound] = findPivot(inData, i, j, member)

pivot = 0;
pivotFound = 0;
a = inData(i).(member);
b = inData(floor((i+j) / 2)).(member);
c = inData(j).(member);

[a, b, c] = o3(a, b, c);

if a < b
    pivot = b;
    pivotFound = 1;
    return;
end

if b < c
    pivot = c;
    pivotFound = 1;
    return;
end

for p = i+1 : j
    if inData(p).(member) ~= inData(i).(member)
        if inData(p).(member) < inData(i).(member)
            pivot = inData(i).(member);
        else
            pivot = inData(p).(member);
        end
        pivotFound = 1;
        return;
    end
    pivotFound = 0;
end
    

%---------------------------

function[a, b, c] = o3(x, y, z)
[x, y] = o2(x, y);
[x, z] = o2(x, z);
[y, z] = o2(y, z);

a = x;
b = y;
c = z;

%---
function[x, y] = o2(x, y)
if (x > y)
    [x, y] = swap(x, y);
end
%---
function [x, y] = swap(x, y)
	t = x;
	x = y;
	y = t;
%----------------------------------
function [inData, p] = partition(inData,i, j, pivot, member)

while i <= j
    while inData(i).(member) < pivot
        i = i + 1;
    end
    while inData(j).(member) >= pivot
        j = j - 1;
    end
    if i < j
        [inData(i), inData(j)] = swap(inData(i), inData(j));
        i = i + 1;
        j = j - 1;
    end
end
p = i;