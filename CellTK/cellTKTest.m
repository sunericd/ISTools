classdef cellTKTest < matlab.unittest.TestCase
    
    properties
    end
    
    methods (Test)
        function testCellCount(testCase)
            % Testing CellCountInitial
            run CellCountInitial.m;
            cellinitialnum = CellCountInitial('Test_Data/cellcount/', '*.tif', 3);
            testCase.assertEqual(cellinitialnum, 8);
            % Testing CellCount
            run CellCount.m;
            cellcountnums = CellCount('Test_Data/cellcount/', '*.tif', 3);
            testCase.assertEqual(cellcountnums(1), 8);
            testCase.assertEqual(cellcountnums(2), 12);
        end
        
        function testCellLite(testCase)
            % Testing CellLite
            run CellLite.m
            cellscore = CellLite('Test_Data/celllite_1', 'Test_Data/celllite_2', '*.jpg', 2)
            testCase.verifyEqual(cellscore, 0.4658)
        end
            
    end
    
end

