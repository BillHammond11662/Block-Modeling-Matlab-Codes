function s=uppercase(S)
% function s=uppercase(S)
%
% moves all chars in string to upper case

if strmatch(class(S),'cell')
 	vals=double(char(S));
	islow = [vals>=97 & vals <=122];
	s = char(vals - islow*32);
    s = cellstr(s);   
else
    vals=double(S);
	islow = [vals>=97 & vals <=122];
	s = char(vals - islow*32);
end;
