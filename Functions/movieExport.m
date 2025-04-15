function [] = movieExport(frames, filename, exportType, dtGif)
    % movieExport
    %
    % Parameters
    % ----------
    % frames : frame
    %   プロットデータ
    % filename : String
    %   ファイル名
    % exportType : String
    %   ファイルの種類 "mp4" or "gif"
    % dtGif : numeric
    %   gifのfps (frame per seconds)
    %   mp4の場合はdtGif=0でよい
    % 
    % Returns
    % -------
    % None
    %

    % 'mp4'のときはdtGif=0でよい
    switch exportType % 今回は出力方法を exporttype 変数で指定することにする
        case 'mp4' % 普通の動画の場合
            video = VideoWriter(filename, 'MPEG-4'); % ファイル名や出力形式などを設定
            open(video);                    % 書き込むファイルを開く
            writeVideo(video, frames);      % ファイルに書き込む
            close(video);                   % 書き込むファイルを閉じる
        case 'gif'
        for i = 1:length(frames)
            [A, map] = rgb2ind(frame2im(frames(i)), 256); % 画像形式変換
            if i == 1
                imwrite(A, map, filename, 'gif', 'DelayTime', 1/dtGif); % 出力形式(FPS)を設定
            else
                imwrite(A, map, filename, 'gif', 'DelayTime', 1/dtGif, 'WriteMode', 'append'); % 2フレーム目以降は"追記"の設定も必要
            end
        end
    end
end