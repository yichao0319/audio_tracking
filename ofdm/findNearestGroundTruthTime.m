function index=findNearestGroundTruthTime(timeRef, groundTruthTime, refIndex)

index=refIndex;
if groundTruthTime(refIndex)>timeRef
    while (groundTruthTime(index)>timeRef)&&(index>1)
        index=index-1;
    end
    if abs(groundTruthTime(index+1)-timeRef)<abs(groundTruthTime(index)-timeRef)
        index=index+1;
    end
elseif groundTruthTime(refIndex)<timeRef
    while (groundTruthTime(index)<timeRef)&&(index<=numel(groundTruthTime))
        index=index+1;
    end
    if abs(groundTruthTime(index-1)-timeRef)<abs(groundTruthTime(index)-timeRef)
        index=index-1;
    end
end