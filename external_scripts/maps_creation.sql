INSERT into Map
SELECT b.binID, g.gid
FROM Bin b
JOIN Gene g ON b.chrom = g.chrom
WHERE (b.binID,g.gid) NOT IN (
   SELECT b.binID, g.gid
   FROM Bin b
   JOIN Gene g ON b.chrom = g.chrom
   WHERE b.end <= g.start OR b.start >= g.end
);
