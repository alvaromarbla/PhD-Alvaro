% Test_CL_Model.m
function Test_CL_Model()
    alpha = deg2rad(5);
    expected_CL = 0.582 + 3.345731*alpha; % First two terms
    assert(abs(CL_Model(alpha) - expected_CL) < 1e-6, 'CL model mismatch');
end