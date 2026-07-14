function str = mat2latex(A, opts)
%MAT2LATEX LaTeX-string van een symbolische of numerieke matrix met [] haken.
%   str = MAT2LATEX(A) zet de matrix A (symbolisch of numeriek) om naar een
%   LaTeX-geformatteerde string waarbij de matrix tussen vierkante haken []
%   staat in plaats van de ronde haken () die latex() standaard gebruikt.
%
%   str = MAT2LATEX(A, math="display") omringt het resultaat met \[ ... \]
%   zodat het rechtstreeks in een LaTeX-document (bv. IguanaTex) geplakt kan
%   worden. Gebruik math="inline" voor $ ... $. Standaard ("none") worden
%   geen wiskunde-scheidingstekens toegevoegd; de string moet dan zelf nog
%   binnen math-mode geplaatst worden.
%
%   str = MAT2LATEX(A, shorthand=true) kort daarnaast cosinussen en sinussen
%   af volgens de gangbare robotica-conventie:
%       cos(th1)     -> c1        sin(th1)     -> s1
%       cos(th2+th3) -> c23       sin(th2+th3) -> s23
%   De labels worden gevormd uit de cijfers in het argument; bevat het
%   argument geen cijfers (bv. cos(alpha)), dan wordt de variabelenaam
%   gebruikt (c_alpha). Niet-herkenbare argumenten blijven ongewijzigd.
%
%   Variabelen genaamd th, th1, th2, ... worden automatisch weergegeven als
%   \theta, \theta_1, \theta_2, ... Andere namen die met "th" beginnen
%   (bv. thickness) blijven ongewijzigd.
%
%   Numerieke matrices worden NIET naar sym omgezet maar rechtstreeks
%   geformatteerd, zodat doubles geen breuken worden (0.1 blijft 0.1 en
%   wordt geen \frac{1}{10}). Het getalformaat is instelbaar via
%   str = MAT2LATEX(A, fmt="%.3f"); standaard is "%.4g". Wetenschappelijke
%   notatie (1.5e+06) wordt omgezet naar 1.5 \cdot 10^{6}. De optie
%   shorthand heeft op numerieke invoer uiteraard geen effect.
%
%   Voorbeeld:
%
%   syms th1 th2 real
%   R = [cos(th1) -sin(th1); sin(th1) cos(th1)];
%   mat2latex(R)
%   mat2latex(R, shorthand=true)   % geeft c_1, s_1 i.p.v. cos(th1), sin(th1)
%   mat2latex(subs(R, th1, 0.3))   % numeriek: 0.9553 i.p.v. een breuk

arguments
    A
    opts.shorthand (1,1) logical = false
    opts.math (1,1) string {mustBeMember(opts.math, ["none", "inline", "display"])} = "none"
    opts.fmt (1,1) string = "%.4g"
end

if isnumeric(A)
    str = numeric2latex(A, char(opts.fmt));
else
    if ~isa(A, 'sym')
        A = sym(A);
    end
    A = renameThetaVars(A);
    if opts.shorthand
        A = applyShorthand(A);
    end
    str = latex(A);
    % latex() zet een matrix in \left(\begin{array}...\end{array}\right);
    % vervang enkel die buitenste array-haken door vierkante haken. Ronde
    % haken binnen elementen (bv. \cos\left(x\right)) blijven ongemoeid.
    str = strrep(str, '\left(\begin{array}', '\left[\begin{array}');
    str = strrep(str, '\end{array}\right)', '\end{array}\right]');
end

switch opts.math
    case "inline"
        str = ['$' str '$'];
    case "display"
        str = ['\[' newline str newline '\]'];
end
end

function str = numeric2latex(A, fmt)
% Formatteer een numerieke matrix rechtstreeks als LaTeX-array, in dezelfde
% stijl als latex(sym(...)) maar zonder sym-omzetting.
elems = arrayfun(@(x) formatNumber(x, fmt), A, 'UniformOutput', false);
if isscalar(elems)
    str = elems{1};
    return
end
rows = cell(size(A, 1), 1);
for r = 1:size(A, 1)
    rows{r} = strjoin(elems(r, :), ' & ');
end
str = ['\left[\begin{array}{' repmat('c', 1, size(A, 2)) '}' newline ...
       strjoin(rows, ['\\\\' newline]) newline ...  % strjoin interpreteert escapes: \\\\ -> \\
       '\end{array}\right]'];
end

function s = formatNumber(x, fmt)
% Eén getal naar LaTeX; Inf/NaN krijgen hun LaTeX-vorm en wetenschappelijke
% notatie wordt omgezet naar \cdot 10^{...}.
if isnan(x)
    s = '\mathrm{NaN}';
    return
elseif isinf(x)
    if x > 0
        s = '\infty';
    else
        s = '-\infty';
    end
    return
end
if isreal(x)
    s = sprintf(fmt, x);
else
    imPart = sprintf(fmt, imag(x));
    if ~startsWith(imPart, '-')
        imPart = ['+' imPart];
    end
    s = [sprintf(fmt, real(x)) imPart '\,\mathrm{i}'];
end
% 1.5e+06 -> 1.5 \cdot 10^{6}, 2e-05 -> 2 \cdot 10^{-5}
s = regexprep(s, 'e\+?(-?)0*(\d+)', ' \\cdot 10^{$1$2}');
end

function A = renameThetaVars(A)
% Hernoem variabelen th, th1, th2, ... naar theta, theta1, theta2, ... zodat
% latex() ze als \theta, \theta_1, ... weergeeft. Andere namen die toevallig
% met "th" beginnen (bv. thickness) blijven ongewijzigd.
for v = reshape(symvar(A), 1, [])
    tok = regexp(char(v), '^th(\d*)$', 'tokens', 'once');
    if ~isempty(tok)
        A = subs(A, v, sym(['theta' tok{1}]));
    end
end
end

function A = applyShorthand(A)
% Vervang cos(arg)/sin(arg) door korte symbolen c<label>/s<label>.
for fn = ["cos" "sin"]
    prefix = extractBefore(fn, 2);          % "c" of "s"
    calls = unique(findSymType(A, fn));
    for k = 1:numel(calls)
        args = children(calls(k));
        arg = args{1};
        numbers = regexp(char(arg), '\d+', 'match');
        if ~isempty(numbers)
            label = strjoin(numbers, '');   % th1 -> 1, th2+th3 -> 23
        elseif isscalar(symvar(arg)) && isequal(arg, symvar(arg))
            label = ['_' char(arg)];        % alpha -> _alpha
        else
            continue                        % geen bruikbare korte naam
        end
        A = subs(A, calls(k), sym(prefix + string(label)));
    end
end
end
