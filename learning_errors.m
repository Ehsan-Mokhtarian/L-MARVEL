function [extra_edges,missing_edges,precision,recall,skeleton_F1_score]=...
    learning_errors(A,B)
    
    A = (A+A')>0;
    B = (B+B')>0;
    skeleton_errors = B-A;
    extra_edges = nnz(skeleton_errors>0)/2;
    missing_edges = nnz(skeleton_errors<0)/2;
    precision = nnz(A.*B)/nnz(B);
    recall = nnz(A.*B)/nnz(A);
    skeleton_F1_score = 2*precision*recall/(precision+recall);
end


