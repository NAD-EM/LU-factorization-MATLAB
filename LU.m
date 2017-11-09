%Calculates solution for a linear system of equations using LU factorization.
%Usage example: "s = LU(magic(5),[1 -1 3 2 3])"

%{
	LU: script to calculate solution for a linear system of equations using LU factorization.
    Copyright (C) 2017  NAD-EM

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    You can contact the author through email:
	<git.nad.em.00@gmail.com>.
	Or through their GitHub profile:
	<https://github.com/NAD-EM>.
%}

function s = LU(A, res)
format long;

s = vec2mat(res, 1);
[m,n] = size(A);
L = eye(m, n);

%matrix A is transformed into matrix U

for column = 1 : n - 1
    for row = column + 1 : m
        Lmn = A(row, column) / A((column + 1) - 1, column);
        L(row, column) = Lmn;
        A(row, :) = A(row, :) + (A((column + 1) - 1, :) * (-Lmn));
    end
end

[s, ~] = linsolve(L, s);
[s, ~] = linsolve(A, s);
end