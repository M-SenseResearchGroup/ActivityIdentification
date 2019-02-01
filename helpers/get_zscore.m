function z = get_zscore(x,mu,sigma)
    z = zeros(size(x));
    for row_ind = 1:size(z,1)
        z(row_ind,:) = (x(row_ind,:)-mu)./sigma; %convert to zscore based on mu and sigma
    end
end