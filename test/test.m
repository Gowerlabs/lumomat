% Execute lumoat test suite under unit-testing framework
%
%
%   (C) Gowerlabs Ltd., 2022
%

% Run low-level API tests
runtests('test_data.m')
runtests('test_SNIRF.m')
runtests('test_lufr.m')
runtests('test_layout_merge.m')


